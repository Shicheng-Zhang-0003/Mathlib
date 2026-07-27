#!/usr/bin/env python3
"""
04_redteam_p0_2_round_correct.py

Run from the folder that CONTAINS the v11S working folder.

Expected layout:

    parent_folder/
        04_redteam_p0_2_round_correct.py
        v11S/
            src/
            tests/
            docs/
            ...

This script applies Red Team P2 item P0-2:

    Replace ml_round() with a correctly-rounded implementation.

It modifies:

    v11S/src/core.c

It adds:

    v11S/tests/test_edge_redteam_p0_2.c

It appends a closure log entry to:

    v11S/docs/CLOSURE_PUNCHLIST.md

Usage:

    python3 04_redteam_p0_2_round_correct.py

Force overwrite:

    python3 04_redteam_p0_2_round_correct.py --force
"""

from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path


NEW_ROUND = r"""ML_API double ml_round(double x) {
    /* MATHLIB_REDP2_P0_2_ROUND_CORRECT */
    if (ml_isnan(x) || ml_isinf(x)) return x;
    if (x == 0.0) return x;

    int neg = ml_signbit(x);
    ml_fp_parts_t p = ml_fp_decompose(neg ? -x : x);

    if (p.kind == ML_FP_ZERO) {
        return x;
    }

    /*
     * If exponent >= 0, the value is already an integer.
     */
    if (p.exp >= 0) {
        return x;
    }

    /*
     * If |x| < 0.5, round to signed zero.
     */
    if (p.exp <= -54) {
        return neg ? -0.0 : 0.0;
    }

    unsigned frac_bits = (unsigned)(-p.exp);

    uint64_t int_part = p.sig >> frac_bits;
    uint64_t frac_mask = (1ULL << frac_bits) - 1ULL;
    uint64_t frac_part = p.sig & frac_mask;
    uint64_t half = 1ULL << (frac_bits - 1);

    /*
     * Round half away from zero.
     *
     * This avoids the broken x +/- 0.5 trick, which can misround
     * values just below halves due to floating-point addition rounding.
     */
    if (frac_part >= half) {
        int_part++;
    }

    double r = (double)int_part;
    return neg ? -r : r;
}
"""


TEST_C = r"""/* MATHLIB_REDP2_P0_2_TEST */
/* Red Team P2 P0-2: ml_round correctness regression tests */
#include "test_harness.h"
#include "ml_core.h"

int main(void) {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Edge Red Team P0-2");

    /*
     * Basic half-away-from-zero behavior.
     */
    ASSERT_TRUE(&ctx, ml_round(0.5) == 1.0, "round(0.5) == 1");
    ASSERT_TRUE(&ctx, ml_round(-0.5) == -1.0, "round(-0.5) == -1");

    ASSERT_TRUE(&ctx, ml_round(1.5) == 2.0, "round(1.5) == 2");
    ASSERT_TRUE(&ctx, ml_round(-1.5) == -2.0, "round(-1.5) == -2");

    ASSERT_TRUE(&ctx, ml_round(2.5) == 3.0, "round(2.5) == 3");
    ASSERT_TRUE(&ctx, ml_round(-2.5) == -3.0, "round(-2.5) == -3");

    /*
     * The old x + 0.5 implementation could misround values just below 0.5.
     *
     * 0x1.fffffffffffffp-2 is the double immediately below 0.5.
     * Correct rounding must return 0, not 1.
     */
    ASSERT_TRUE(&ctx,
        ml_round(0x1.fffffffffffffp-2) == 0.0,
        "round(nextbelow(0.5)) == 0");

    double rn = ml_round(-0x1.fffffffffffffp-2);
    ASSERT_TRUE(&ctx,
        rn == 0.0 && ml_signbit(rn) != 0,
        "round(-nextbelow(0.5)) == -0");

    /*
     * Values just above 0.5 must round to 1 / -1.
     */
    ASSERT_TRUE(&ctx,
        ml_round(0x1.0000000000001p-1) == 1.0,
        "round(nextabove(0.5)) == 1");

    ASSERT_TRUE(&ctx,
        ml_round(-0x1.0000000000001p-1) == -1.0,
        "round(-nextabove(0.5)) == -1");

    /*
     * Signed zero behavior.
     */
    ASSERT_TRUE(&ctx,
        ml_round(0.0) == 0.0 && ml_signbit(ml_round(0.0)) == 0,
        "round(+0) preserves +0");

    ASSERT_TRUE(&ctx,
        ml_round(-0.0) == 0.0 && ml_signbit(ml_round(-0.0)) != 0,
        "round(-0) preserves -0");

    /*
     * Non-finite behavior.
     */
    ASSERT_TRUE(&ctx,
        ml_isnan(ml_round(ml_make_nan())),
        "round(NaN) is NaN");

    ASSERT_TRUE(&ctx,
        ml_isinf(ml_round(ml_make_inf(0))),
        "round(+inf) is +inf");

    ASSERT_TRUE(&ctx,
        ml_isinf(ml_round(-ml_make_inf(0))),
        "round(-inf) is -inf");

    /*
     * Large integers remain unchanged.
     */
    ASSERT_TRUE(&ctx,
        ml_round(1e16) == 1e16,
        "round(large integer) identity");

    ASSERT_TRUE(&ctx,
        ml_round(-1e16) == -1e16,
        "round(-large integer) identity");

    /*
     * Boundary near 2^52.
     *
     * 0x1.fffffffffffffp+51 == 2^52 - 0.5
     * It must round to 2^52.
     */
    ASSERT_TRUE(&ctx,
        ml_round(0x1.fffffffffffffp+51) == 0x1p+52,
        "round(2^52 - 0.5) == 2^52");

    ASSERT_TRUE(&ctx,
        ml_round(-0x1.fffffffffffffp+51) == -0x1p+52,
        "round(-(2^52 - 0.5)) == -2^52");

    return ml_test_summary(&ctx);
}
"""


PUNCHLIST_LOG = r"""
<!-- MATHLIB_REDP2_P0_2_LOG -->
## Red Team P2 P0-2

- Replaced `ml_round()` with exact decomposition-based round-half-away-from-zero.
- Fixed misrounding near half-integers caused by the old `x +/- 0.5` trick.
- Added regression suite `tests/test_edge_redteam_p0_2.c`.
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
    if (root / "src" / "core.c").is_file() and (root / "include" / "mathlib").is_dir():
        print("[note] Running from inside v11S; treating current directory as v11S.")
        return root.parent, root

    fail(
        "Run this script from the folder that CONTAINS the v11S directory, "
        "or from inside v11S itself."
    )


def patch_core(v11s: Path, force: bool) -> None:
    path = v11s / "src" / "core.c"

    if not path.is_file():
        fail(f"Missing expected file: {path}")

    text = normalize(path.read_text(encoding="utf-8"))

    marker = "MATHLIB_REDP2_P0_2_ROUND_CORRECT"

    if marker in text and not force:
        print(f"[skip] {path}: P0-2 marker already present")
        return

    pattern = re.compile(
        r"(?ms)^[ \t]*ML_API[ \t]+double[ \t]+ml_round\("
        r"double[ \t]+x\)[ \t]*\{"
        r".*?"
        r"(?=^[ \t]*#ifndef[ \t]+ML_PI|\Z)"
    )

    patched, count = pattern.subn(
        lambda m: NEW_ROUND + "\n",
        text,
        count=1,
    )

    if count != 1:
        fail(
            f"{path}: expected exactly one ml_round() definition match, "
            f"got {count}. Source may have drifted."
        )

    write_text(path, patched)


def write_test(v11s: Path, force: bool) -> None:
    path = v11s / "tests" / "test_edge_redteam_p0_2.c"

    if path.exists() and not force:
        try:
            old = normalize(path.read_text(encoding="utf-8"))
            if "MATHLIB_REDP2_P0_2_TEST" in old:
                print(f"[skip] {path}: already present")
                return
        except UnicodeDecodeError:
            pass

    write_text(path, TEST_C)


def append_log(v11s: Path) -> None:
    path = v11s / "docs" / "CLOSURE_PUNCHLIST.md"

    marker = "<!-- MATHLIB_REDP2_P0_2_LOG -->"

    if not path.exists():
        write_text(
            path,
            "# MathLib v11S Closure Punchlist\n" + PUNCHLIST_LOG,
        )
        return

    text = normalize(path.read_text(encoding="utf-8"))

    if marker in text:
        print(f"[skip] {path}: P0-2 log already present")
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
    print("  Red Team P2 P0-2: ml_round correct rounding fix")
    print("=========================================================")
    print(f"Root: {root}")
    print(f"v11S: {v11s}")
    print(f"Force: {force}")
    print("---------------------------------------------------------")

    patch_core(v11s, force)
    write_test(v11s, force)
    append_log(v11s)
    archive_self(v11s, force)

    print("---------------------------------------------------------")
    print("P0-2 changes applied.")
    print("")
    print("Next verification step:")
    print("")
    print("    cd v11S")
    print("    MATHLIB_EDGE_SANITIZERS=1 bash tests/run_edge_tests.sh")
    print("")
    print("Expected new suite:")
    print("")
    print("    --- test_edge_redteam_p0_2 ---")
    print("    [Edge Red Team P0-2] Passed: ..., Failed: 0")
    print("")
    print("After that passes, the next block is:")
    print("")
    print("    P0-3: fix ml_exp() overflow/underflow thresholds")
    print("=========================================================")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
