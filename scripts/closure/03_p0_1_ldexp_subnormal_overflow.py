#!/usr/bin/env python3
"""
03_p0_1_ldexp_subnormal_overflow.py

Run from the folder that CONTAINS the v11S working folder.

Expected layout:

    parent_folder/
        03_p0_1_ldexp_subnormal_overflow.py
        v11S/
            src/
            tests/
            docs/
            ...

This script implements Red Team P2 item P0-1:

    Fix ml_ldexp_pure() false overflow for subnormal significands.

It modifies:

    v11S/src/core.c

It adds:

    v11S/tests/test_edge_redteam_p0_1.c

It appends a closure log entry to:

    v11S/docs/CLOSURE_PUNCHLIST.md

Usage:

    python3 03_p0_1_ldexp_subnormal_overflow.py

Force overwrite:

    python3 03_p0_1_ldexp_subnormal_overflow.py --force
"""

from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path


NEW_LDEXP = r'''ML_API double ml_ldexp_pure(double x, int exp) {
    /* MATHLIB_CLOSURE_P2_P0_1_LDEXP_NORMALIZE */
    ml_fp_parts_t p = ml_fp_decompose(x);

    if (p.kind == ML_FP_ZERO || p.kind == ML_FP_INF || p.kind == ML_FP_NAN) {
        return x;
    }

    long long e = (long long)p.exp + (long long)exp;
    uint64_t sig = p.sig;

    if (sig == 0) {
        return ml_copysign(0.0, x);
    }

    /*
     * Normalize subnormal significands BEFORE applying overflow bounds.
     *
     * The old code applied the normal-number overflow bound (e > 971)
     * directly to subnormal significands, causing false overflow for
     * small sig + large positive exponent.
     *
     * Example:
     *   smallest subnormal = 2^-1074
     *   ldexp(smallest subnormal, 2046) = 2^972, which is finite.
     *
     * The old path saw:
     *   new_exp = -1074 + 2046 = 972
     *   972 > 971 -> infinity
     *
     * That was incorrect.
     */
    while (sig < (1ULL << 52) && e > -2098) {
        sig <<= 1;
        e--;
    }

    if (e > 971) {
        return ml_make_inf(p.sign);
    }

    if (e < -2098) {
        uint64_t bits = p.sign ? 0x8000000000000000ULL : 0ULL;
        double d;
        memcpy(&d, &bits, sizeof(double));
        return d;
    }

    return ml_fp_compose(sig, (int)e, p.sign);
}
'''


TEST_C = r'''/* MATHLIB_CLOSURE_P2_P0_1_TEST */
/* Red Team P2 P0-1: ldexp subnormal overflow regression tests */
#include "test_harness.h"
#include "ml_core.h"

#include <stdint.h>
#include <string.h>

int main(void) {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Edge Red Team P0-1");

    uint64_t min_bits = 1;
    double min_sub;
    memcpy(&min_sub, &min_bits, sizeof(double));

    /*
     * Subnormal + large positive exponent must not falsely overflow.
     */
    double r2046 = ml_ldexp_pure(min_sub, 2046);
    ASSERT_TRUE(&ctx,
        ml_isfinite(r2046),
        "ldexp(min_subnormal, 2046) is finite");

    ASSERT_TRUE(&ctx,
        r2046 == 0x1p972,
        "ldexp(min_subnormal, 2046) == 2^972");

    double rneg = ml_ldexp_pure(-min_sub, 2046);
    ASSERT_TRUE(&ctx,
        ml_isfinite(rneg) && ml_signbit(rneg) != 0,
        "ldexp(-min_subnormal, 2046) preserves sign");

    ASSERT_TRUE(&ctx,
        rneg == -0x1p972,
        "ldexp(-min_subnormal, 2046) == -2^972");

    /*
     * Finite boundary near the top of the double range.
     */
    double r2097 = ml_ldexp_pure(min_sub, 2097);
    ASSERT_TRUE(&ctx,
        ml_isfinite(r2097),
        "ldexp(min_subnormal, 2097) is finite");

    ASSERT_TRUE(&ctx,
        r2097 == 0x1p1023,
        "ldexp(min_subnormal, 2097) == 2^1023");

    double r2098 = ml_ldexp_pure(min_sub, 2098);
    ASSERT_TRUE(&ctx,
        ml_isinf(r2098) && r2098 > 0.0,
        "ldexp(min_subnormal, 2098) overflows to +inf");

    /*
     * Normal-input sanity boundaries.
     */
    double n1023 = ml_ldexp_pure(1.0, 1023);
    ASSERT_TRUE(&ctx,
        ml_isfinite(n1023),
        "ldexp(1, 1023) is finite");

    ASSERT_TRUE(&ctx,
        n1023 == 0x1p1023,
        "ldexp(1, 1023) == 2^1023");

    ASSERT_TRUE(&ctx,
        ml_isinf(ml_ldexp_pure(1.0, 1024)) && ml_ldexp_pure(1.0, 1024) > 0.0,
        "ldexp(1, 1024) overflows to +inf");

    /*
     * Existing underflow behavior must remain intact.
     */
    ASSERT_TRUE(&ctx,
        ml_ldexp_pure(1.0, -1074) == min_sub,
        "ldexp(1, -1074) == min_subnormal");

    ASSERT_TRUE(&ctx,
        ml_ldexp_pure(1.0, -1075) == 0.0,
        "ldexp(1, -1075) underflows to zero");

    ASSERT_TRUE(&ctx,
        ml_signbit(ml_ldexp_pure(-1.0, -1075)) != 0,
        "ldexp(-1, -1075) preserves negative zero sign");

    return ml_test_summary(&ctx);
}
'''


LOG_ENTRY = r'''
<!-- MATHLIB_CLOSURE_P2_P0_1_LOG -->
## Red Team P2 P0-1

- Fixed `ml_ldexp_pure()` false overflow for subnormal significands.
- Normalization now occurs before overflow judgment.
- Added regression suite `tests/test_edge_redteam_p0_1.c`.
'''


def fail(message: str) -> None:
    print("ERROR: " + message)
    sys.exit(1)


def write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(content)
    print(f"[write] {path}")


def normalize(text: str) -> str:
    return text.replace("\r\n", "\n").replace("\r", "\n")


def patch_core(v11s: Path, force: bool) -> None:
    core_path = v11s / "src" / "core.c"

    if not core_path.is_file():
        fail(f"Missing expected file: {core_path}")

    original = core_path.read_text(encoding="utf-8")
    text = normalize(original)

    marker = "MATHLIB_CLOSURE_P2_P0_1_LDEXP_NORMALIZE"

    if marker in text and not force:
        print(f"[skip] {core_path}: P0-1 marker already present")
        return

    pattern = re.compile(
        r"(?ms)^ML_API double ml_ldexp_pure\(double x, int exp\)\s*\{.*?^\}[ \t]*(?:\n|\Z)"
    )

    patched, count = pattern.subn(lambda m: NEW_LDEXP, text, count=1)

    if count != 1:
        fail(
            f"{core_path}: expected exactly one ml_ldexp_pure() match, got {count}. "
            "Source may have drifted."
        )

    write_text(core_path, patched)


def write_test(v11s: Path, force: bool) -> None:
    test_path = v11s / "tests" / "test_edge_redteam_p0_1.c"

    if test_path.exists() and not force:
        try:
            old = test_path.read_text(encoding="utf-8")
            if "MATHLIB_CLOSURE_P2_P0_1_TEST" in old:
                print(f"[skip] {test_path}: already present")
                return
        except UnicodeDecodeError:
            pass

    write_text(test_path, TEST_C)


def append_log(v11s: Path) -> None:
    punch_path = v11s / "docs" / "CLOSURE_PUNCHLIST.md"

    if not punch_path.exists():
        print(f"[skip] {punch_path}: not found")
        return

    text = normalize(punch_path.read_text(encoding="utf-8"))

    if "<!-- MATHLIB_CLOSURE_P2_P0_1_LOG -->" in text:
        print(f"[skip] {punch_path}: P0-1 log already present")
        return

    with open(punch_path, "a", encoding="utf-8", newline="\n") as fh:
        fh.write(LOG_ENTRY)

    print(f"[append] {punch_path}")


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

    root = Path.cwd()
    v11s = root / "v11S"

    if not v11s.is_dir():
        fail("Run this script from the folder that CONTAINS the v11S directory.")

    print("=========================================================")
    print("  Red Team P2 P0-1: ml_ldexp_pure subnormal overflow fix")
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
    print("P0-1 changes applied.")
    print("")
    print("Next verification step:")
    print("")
    print("    cd v11S")
    print("    MATHLIB_EDGE_SANITIZERS=1 bash tests/run_edge_tests.sh")
    print("")
    print("The new test will be picked up automatically because it matches:")
    print("")
    print("    tests/test_edge_*.c")
    print("=========================================================")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
