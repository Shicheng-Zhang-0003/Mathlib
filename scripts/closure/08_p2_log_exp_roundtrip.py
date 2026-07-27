#!/usr/bin/env python3
"""
08_p2_log_exp_roundtrip.py

Run from the folder that CONTAINS the v11S working folder.

Expected layout:

    parent_folder/
        08_p2_log_exp_roundtrip.py
        v11S/
            src/
            tests/
            docs/
            ...

This script fixes the Red Team P2 issue:

    fuzz_god_mode log(exp(x)) roundtrip tolerance is too strict.

The old fuzzer used:

    ml_fabs(log_exp - x) < 1e-13

But for x in [512, 1024), one double ULP is about 1.14e-13,
so the old absolute tolerance is tighter than one ULP.

This script:

1. Adds a ULP-distance helper to fuzz_god_mode.c
2. Replaces the log(exp(x)) tolerance with a ULP-aware tolerance
3. Improves ml_log() final reconstruction using ML_FMA
4. Adds a deterministic regression edge test
5. Appends a closure log entry

Targets:

    v11S/tests/fuzz_god_mode.c
    v11S/src/exp_log.c
    v11S/tests/test_edge_redteam_p2_logexp.c
    v11S/docs/CLOSURE_PUNCHLIST.md

Usage:

    python3 08_p2_log_exp_roundtrip.py

Force overwrite:

    python3 08_p2_log_exp_roundtrip.py --force
"""

from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path


ULP_HELPER = r"""/* MATHLIB_CLOSURE_P2_LOG_EXP_ROUNDTRIP_HELPER */
static uint64_t ml_fuzz_ulp_distance(double a, double b) {
    uint64_t ia, ib;

    memcpy(&ia, &a, sizeof(uint64_t));
    memcpy(&ib, &b, sizeof(uint64_t));

    if (ia >> 63) ia = 0x8000000000000000ULL - ia;
    if (ib >> 63) ib = 0x8000000000000000ULL - ib;

    return ia > ib ? ia - ib : ib - ia;
}
"""


NEW_LOGEXP_BLOCK = r"""if (x < 700.0) {
    double exp_x = ml_exp(x);
    double log_exp = ml_log(exp_x);
    double diff = ml_fabs(log_exp - x);
    uint64_t ulps = ml_fuzz_ulp_distance(log_exp, x);

    /*
     * MATHLIB_CLOSURE_P2_LOG_EXP_ROUNDTRIP_TOLERANCE
     *
     * The old absolute tolerance 1e-13 is smaller than one double ULP
     * for x in [512, 1024), where one ULP is approximately 1.14e-13.
     *
     * Therefore a roundtrip error of one ULP could fail the old test
     * even though the result is numerically sane.
     *
     * Accept either:
     *   - very small absolute error, or
     *   - <= 8 ULP roundtrip error.
     */
    if (ml_isnan(log_exp) || ml_isinf(log_exp) || !(diff <= 1e-13 || ulps <= 8)) {
        failed++;
        printf("FAIL: log(exp(x)) | x: %.17g | exp(x): %.17g | log(exp(x)): %.17g | diff: %.17g | ulps: %llu\n",
               x, exp_x, log_exp, diff, (unsigned long long)ulps);
    } else {
        passed++;
    }
}
"""


EDGE_TEST_C = r"""/* MATHLIB_CLOSURE_P2_LOG_EXP_ROUNDTRIP_TEST */
/* Red Team P2: log(exp(x)) roundtrip regression tests */
#include "test_harness.h"
#include "ml_core.h"
#include "ml_exp_log.h"

#include <stdint.h>
#include <string.h>

static uint64_t redteam_ulp_distance(double a, double b) {
    uint64_t ia, ib;

    memcpy(&ia, &a, sizeof(uint64_t));
    memcpy(&ib, &b, sizeof(uint64_t));

    if (ia >> 63) ia = 0x8000000000000000ULL - ia;
    if (ib >> 63) ib = 0x8000000000000000ULL - ib;

    return ia > ib ? ia - ib : ib - ia;
}

static void check_roundtrip(ml_test_ctx_t *ctx, double x) {
    double e = ml_exp(x);
    double l = ml_log(e);

    ASSERT_TRUE(ctx,
        ml_isfinite(e) && ml_isfinite(l),
        "exp/log roundtrip finite");

    if (!ml_isfinite(e) || !ml_isfinite(l)) {
        return;
    }

    double diff = ml_fabs(l - x);
    uint64_t ulps = redteam_ulp_distance(l, x);

    ASSERT_TRUE(ctx,
        diff <= 1e-13 || ulps <= 8,
        "log(exp(x)) roundtrip within 1e-13 absolute or 8 ULP");
}

int main(void) {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Edge Red Team P2 Log/Exp Roundtrip");

    /*
     * These values came from deterministic fuzz_god_mode seed 123456789.
     * They are in the double range where 1 ULP is about 1.14e-13.
     */
    check_roundtrip(&ctx, 586.978120536999);
    check_roundtrip(&ctx, 627.643520770428);
    check_roundtrip(&ctx, 671.779855933938);
    check_roundtrip(&ctx, 674.258560721883);
    check_roundtrip(&ctx, 628.762675276381);
    check_roundtrip(&ctx, 588.572237449033);
    check_roundtrip(&ctx, 609.755943347586);
    check_roundtrip(&ctx, 605.736818912316);
    check_roundtrip(&ctx, 688.953238394555);
    check_roundtrip(&ctx, 628.651212727954);
    check_roundtrip(&ctx, 651.39389487514);
    check_roundtrip(&ctx, 672.896168042392);
    check_roundtrip(&ctx, 630.227606571385);
    check_roundtrip(&ctx, 632.544160650434);
    check_roundtrip(&ctx, 654.068812567841);
    check_roundtrip(&ctx, 628.888303241175);
    check_roundtrip(&ctx, 544.08798345555);

    return ml_test_summary(&ctx);
}
"""


PUNCHLIST_LOG = r"""
<!-- MATHLIB_CLOSURE_P2_LOG_EXP_ROUNDTRIP_LOG -->
## Red Team P2 Log/Exp Roundtrip Fix

- Replaced over-tight `log(exp(x))` fuzzer tolerance with ULP-aware tolerance.
- Added `ml_fuzz_ulp_distance()` to `fuzz_god_mode.c`.
- Improved `ml_log()` final reconstruction using `ML_FMA`.
- Added deterministic regression suite `tests/test_edge_redteam_p2_logexp.c`.
"""


def fail(message: str) -> None:
    print("ERROR: " + message)
    sys.exit(1)


def normalize(text: str) -> str:
    return text.replace("\r\n", "\n").replace("\r", "\n")


def write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(content)
    print(f"[write] {path}")


def locate_v11s() -> tuple[Path, Path]:
    root = Path.cwd()

    candidate = root / "v11S"
    if candidate.is_dir():
        return root, candidate

    # Convenience: allow running from inside v11S itself.
    if (root / "src" / "exp_log.c").is_file() and (root / "tests" / "fuzz_god_mode.c").is_file():
        print("[note] Running from inside v11S; treating current directory as v11S.")
        return root.parent, root

    fail(
        "Run this script from the folder that CONTAINS the v11S directory, "
        "or from inside v11S itself."
    )


def patch_fuzzer(v11s: Path, force: bool) -> None:
    path = v11s / "tests" / "fuzz_god_mode.c"

    if not path.is_file():
        fail(f"Missing expected file: {path}")

    text = normalize(path.read_text(encoding="utf-8"))
    changed = False

    helper_marker = "MATHLIB_CLOSURE_P2_LOG_EXP_ROUNDTRIP_HELPER"
    block_marker = "MATHLIB_CLOSURE_P2_LOG_EXP_ROUNDTRIP_TOLERANCE"

    if helper_marker in text:
        print(f"[skip] {path}: ULP helper already present")
    else:
        helper_pattern = re.compile(
            r"(?m)^[ \t]*void test_ieee754_specials\(\)\s*\{"
        )

        patched, count = helper_pattern.subn(
            lambda m: ULP_HELPER + "\n" + m.group(0),
            text,
            count=1,
        )

        if count != 1:
            fail(
                f"{path}: could not find test_ieee754_specials() anchor "
                "for ULP helper insertion."
            )

        text = patched
        changed = True
        print(f"[patch] {path}: inserted ULP helper")

    if block_marker in text and not force:
        print(f"[skip] {path}: log(exp(x)) tolerance block already present")
    else:
        block_pattern = re.compile(
            r"(?s)if \(x < 700\.0\) \{.*?\}\s*(?=double prod = x \* y;)"
        )

        patched, count = block_pattern.subn(
            lambda m: NEW_LOGEXP_BLOCK + "\n",
            text,
            count=1,
        )

        if count != 1:
            fail(
                f"{path}: could not find the log(exp(x)) test block. "
                "Source may have drifted."
            )

        text = patched
        changed = True
        print(f"[patch] {path}: replaced log(exp(x)) tolerance block")

    if changed:
        write_text(path, text)
    else:
        print(f"[skip] {path}: already patched")


def patch_ml_log(v11s: Path) -> None:
    path = v11s / "src" / "exp_log.c"

    if not path.is_file():
        fail(f"Missing expected file: {path}")

    text = normalize(path.read_text(encoding="utf-8"))

    marker = "MATHLIB_CLOSURE_P2_LOG_FMA_RECONSTRUCT"

    if marker in text:
        print(f"[skip] {path}: ml_log FMA reconstruction already present")
        return

    pattern = re.compile(
        r"(?m)^([ \t]*)return z \* poly \+ e \* ML_LN2;[ \t]*$"
    )

    def repl(m: re.Match[str]) -> str:
        indent = m.group(1)
        return (
            f"{indent}/* MATHLIB_CLOSURE_P2_LOG_FMA_RECONSTRUCT */\n"
            f"{indent}return ML_FMA((double)e, ML_LN2, z * poly);"
        )

    patched, count = pattern.subn(repl, text, count=1)

    if count != 1:
        fail(
            f"{path}: could not find ml_log return expression. "
            "Source may have drifted."
        )

    write_text(path, patched)


def write_edge_test(v11s: Path, force: bool) -> None:
    path = v11s / "tests" / "test_edge_redteam_p2_logexp.c"

    if path.exists() and not force:
        try:
            old = normalize(path.read_text(encoding="utf-8"))
            if "MATHLIB_CLOSURE_P2_LOG_EXP_ROUNDTRIP_TEST" in old:
                print(f"[skip] {path}: already present")
                return
        except UnicodeDecodeError:
            pass

    write_text(path, EDGE_TEST_C)


def append_log(v11s: Path) -> None:
    path = v11s / "docs" / "CLOSURE_PUNCHLIST.md"

    if not path.exists():
        write_text(
            path,
            "# MathLib v11S Closure Punchlist\n" + PUNCHLIST_LOG,
        )
        return

    text = normalize(path.read_text(encoding="utf-8"))

    if "<!-- MATHLIB_CLOSURE_P2_LOG_EXP_ROUNDTRIP_LOG -->" in text:
        print(f"[skip] {path}: log/exp roundtrip log already present")
        return

    with open(path, "a", encoding="utf-8", newline="\n") as fh:
        fh.write(PUNCHLIST_LOG)

    print(f"[append] {path}")


def archive_self(v11s: Path, force: bool) -> None:
    try:
        source_script = Path(__file__).resolve()
        archived_script = v11s / "scripts" / "closure" / source_script.name

        if source_script == archived_script:
            return

        if archived_script.exists() and not force:
            print(f"[skip] {archived_script}: already archived")
            return

        archived_script.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source_script, archived_script)
        print(f"[archive] {archived_script}")

    except NameError:
        print("[note] Could not archive script because __file__ is unavailable.")


def main() -> int:
    force = "--force" in sys.argv[1:]

    root, v11s = locate_v11s()

    print("=========================================================")
    print("  Red Team P2: log/exp roundtrip tolerance fix")
    print("=========================================================")
    print(f"Root: {root}")
    print(f"v11S: {v11s}")
    print(f"Force: {force}")
    print("---------------------------------------------------------")

    patch_fuzzer(v11s, force)
    patch_ml_log(v11s)
    write_edge_test(v11s, force)
    append_log(v11s)
    archive_self(v11s, force)

    print("---------------------------------------------------------")
    print("P2 log/exp roundtrip changes applied.")
    print("")
    print("Now rebuild the fuzzer and tests.")
    print("")
    print("If using CMake:")
    print("")
    print("    cd v11S")
    print("    cmake --build build")
    print("    ./build/fuzz_god_mode 123456789")
    print("")
    print("If using Make:")
    print("")
    print("    cd v11S")
    print("    make fuzz_god_mode")
    print("    ./build/fuzz_god_mode 123456789")
    print("")
    print("Also run edge tests:")
    print("")
    print("    MATHLIB_EDGE_SANITIZERS=1 bash tests/run_edge_tests.sh")
    print("=========================================================")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
