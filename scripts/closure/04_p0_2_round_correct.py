#!/usr/bin/env python3
"""
04_p0_2_round_correct.py

Run from the folder that CONTAINS the v11S working folder.

Expected layout:

    parent_folder/
        04_p0_2_round_correct.py
        v11S/
            src/
            tests/
            docs/
            ...

This script implements Red Team P2 item P0-2:

    Replace ml_round() with a correctly-rounded implementation.

It modifies:

    v11S/src/core.c

It adds:

    v11S/tests/test_edge_redteam_p0_2.c

It appends a closure log entry to:

    v11S/docs/CLOSURE_PUNCHLIST.md

Usage:

    python3 04_p0_2_round_correct.py

Force overwrite:

    python3 04_p0_2_round_correct.py --force
"""

from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path


NEW_ROUND = r'''ML_API double ml_round(double x) {
    /* MATHLIB_CLOSURE_P2_P0_2_ROUND_EXACT */
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
'''


TEST_C = r'''/* MATHLIB_CLOSURE_P2_P0_2_TEST */
/* Red Team P2 P0-2: ml_round correctness regression tests */
#include "test_harness.h"
#include "ml_core.h"

int main(void) {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Edge Red Team P0-2");

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
     * Signed zero behavior.
     */
    ASSERT_TRUE(&ctx,
        ml_round(0.0) == 0.0 && ml_signbit(ml_round(0.0)) == 0,
        "round(+0) preserves +0");

    ASSERT_TRUE(&ctx,
        ml_round(-0.0) == 0.0 && ml_signbit(ml_round(-0.0)) != 0,
        "round(-0) preserves -0");

    /*
     * Exact half cases: round half away from zero.
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
     * Large integers remain unchanged.
     */
    ASSERT_TRUE(&ctx,
        ml_round(1e16) == 1e16,
        "round(large integer) identity");

    ASSERT_TRUE(&ctx,
        ml_round(-1e16) == -1e16,
        "round(-large integer) identity");

    return ml_test_summary(&ctx);
}
'''


PUNCHLIST_LOG = r'''
<!-- MATHLIB_CLOSURE_P2_P0_2_LOG -->
## Red Team P2 P0-2

- Replaced `ml_round()` with exact decomposition-based round-half-away-from-zero.
- Fixed misrounding near half-integers caused by the old `x +/- 0.5` trick.
- Added regression suite `tests/test_edge_redteam_p0_2.c`.
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

    marker = "MATHLIB_CLOSURE_P2_P0_2_ROUND_EXACT"

    if marker in text and not force:
        print(f"[skip] {core_path}: P0-2 marker already present")
        return

    pattern = re.compile(
        r"(?ms)^ML_API double ml_round\(double x\)\s*\{.*?^\}[ \t]*(?:\n|\Z)"
    )

    patched, count = pattern.subn(lambda m: NEW_ROUND, text, count=1)

    if count != 1:
        fail(
            f"{core_path}: expected exactly one ml_round() match, got {count}. "
            "Source may have drifted."
        )

    write_text(core_path, patched)


def write_test(v11s: Path, force: bool) -> None:
    test_path = v11s / "tests" / "test_edge_redteam_p0_2.c"

    if test_path.exists() and not force:
        try:
            old = test_path.read_text(encoding="utf-8")
            if "MATHLIB_CLOSURE_P2_P0_2_TEST" in old:
                print(f"[skip] {test_path}: already present")
                return
        except UnicodeDecodeError:
            pass

    write_text(test_path, TEST_C)


def append_log(v11s: Path) -> None:
    punch_path = v11s / "docs" / "CLOSURE_PUNCHLIST.md"

    if not punch_path.exists():
        write_text(
            punch_path,
            "# MathLib v11S Closure Punchlist\n" + PUNCHLIST_LOG,
        )
        return

    text = normalize(punch_path.read_text(encoding="utf-8"))

    if "<!-- MATHLIB_CLOSURE_P2_P0_2_LOG -->" in text:
        print(f"[skip] {punch_path}: P0-2 log already present")
        return

    with open(punch_path, "a", encoding="utf-8", newline="\n") as fh:
        fh.write(PUNCHLIST_LOG)

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
    print("  Red Team P2 P0-2: correct ml_round() implementation")
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
    print("After that passes, the next block is:")
    print("")
    print("    P0-3: fix ml_exp() overflow/underflow thresholds")
    print("=========================================================")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
