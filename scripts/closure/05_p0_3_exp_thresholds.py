#!/usr/bin/env python3
"""
05_p0_3_exp_thresholds.py

Run from the folder that CONTAINS the v11S working folder.

Expected layout:

    parent_folder/
        05_p0_3_exp_thresholds.py
        v11S/
            src/
            tests/
            docs/
            ...

This script implements Red Team P2 item P0-3:

    Fix ml_exp() premature overflow / underflow thresholds.

It modifies:

    v11S/src/exp_log.c

It adds:

    v11S/tests/test_edge_redteam_p0_3.c

It appends a closure log entry to:

    v11S/docs/CLOSURE_PUNCHLIST.md

Usage:

    python3 05_p0_3_exp_thresholds.py

Force overwrite:

    python3 05_p0_3_exp_thresholds.py --force
"""

from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path


MACRO_BLOCK = r"""/* MATHLIB_CLOSURE_P2_P0_3_EXP_LIMITS */
#ifndef ML_LOG_DBL_MAX
#define ML_LOG_DBL_MAX 709.782712893384
#endif

#ifndef ML_LOG_UNDERFLOW
#define ML_LOG_UNDERFLOW (-745.133219101941)
#endif
"""


NEW_EXP = r"""ML_API double ml_exp(double x) {
    /* MATHLIB_CLOSURE_P2_P0_3_EXP_LIMITS */
    if (ml_isnan(x)) return x;
    if (ml_isinf(x)) return (x > 0.0) ? ml_make_inf(0) : 0.0;
    if (x == 0.0) return 1.0;

    /*
     * Use tighter, double-limit-aware thresholds.
     *
     * Old code used:
     *   x > 709.78      -> inf
     *   x < -745.13     -> 0
     *
     * Those thresholds were too conservative and could prematurely
     * overflow or underflow values that are still representable.
     */
    if (x > ML_LOG_DBL_MAX) {
        return ml_make_inf(0);
    }

    if (x < ML_LOG_UNDERFLOW) {
        return 0.0;
    }

    double n = ml_round(x / ML_LN2);
    double r = x - n * 0.69314718036912381649 - n * 1.90821490974462528503e-10;

    static const double inv_fact[] = {
        1.0, 1.0, 0.5, 0.16666666666666666, 0.041666666666666664,
        0.008333333333333333, 0.001388888888888889, 0.0001984126984126984,
        2.48015873015873e-05, 2.7557319223985893e-06, 2.7557319223985888e-07,
        2.505210838544172e-08, 2.08767569878681e-09, 1.6059043836821613e-10,
        1.1470745597729725e-11, 7.647163731819816e-13, 4.779477332387385e-14,
        2.8114572543455206e-15, 1.5619206968586226e-16, 8.22063524662433e-18
    };

    double result = inv_fact[19];
    for (int i = 18; i >= 1; i--) {
        result = ML_FMA(result, r, inv_fact[i]);
    }
    result = ML_FMA(result, r, 1.0);

    return ml_ldexp_pure(result, (int)n);
}
"""


TEST_C = r"""/* MATHLIB_CLOSURE_P2_P0_3_TEST */
/* Red Team P2 P0-3: ml_exp threshold regression tests */
#include "test_harness.h"
#include "ml_core.h"
#include "ml_exp_log.h"

int main(void) {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Edge Red Team P0-3");

    /*
     * Old overflow threshold was too low:
     *
     *   if (x > 709.78) return inf;
     *
     * The true approximate overflow boundary is:
     *
     *   log(DBL_MAX) ~= 709.782712893384
     *
     * Therefore values slightly above 709.78 but below log(DBL_MAX)
     * must remain finite.
     */
    ASSERT_TRUE(&ctx,
        ml_isfinite(ml_exp(709.781)),
        "exp(709.781) must remain finite");

    ASSERT_TRUE(&ctx,
        ml_isfinite(ml_exp(709.782)),
        "exp(709.782) must remain finite");

    /*
     * Values above log(DBL_MAX) must overflow.
     */
    ASSERT_TRUE(&ctx,
        ml_isinf(ml_exp(709.783)) && ml_exp(709.783) > 0.0,
        "exp(709.783) must overflow to +inf");

    ASSERT_TRUE(&ctx,
        ml_isinf(ml_exp(710.0)) && ml_exp(710.0) > 0.0,
        "exp(710) must overflow to +inf");

    /*
     * Old underflow threshold was too high:
     *
     *   if (x < -745.13) return 0;
     *
     * The approximate underflow-to-zero boundary is closer to:
     *
     *   -745.133219101941
     *
     * Therefore -745.132 should still produce a positive subnormal
     * or normal result, not zero.
     */
    ASSERT_TRUE(&ctx,
        ml_exp(-745.132) > 0.0,
        "exp(-745.132) must not underflow prematurely");

    ASSERT_TRUE(&ctx,
        ml_exp(-745.2) == 0.0,
        "exp(-745.2) underflows to zero");

    /*
     * Basic special-case sanity must remain intact.
     */
    ASSERT_TRUE(&ctx,
        ml_exp(0.0) == 1.0,
        "exp(0) == 1");

    ASSERT_TRUE(&ctx,
        ml_isnan(ml_exp(ml_make_nan())),
        "exp(NaN) is NaN");

    ASSERT_TRUE(&ctx,
        ml_isinf(ml_exp(ml_make_inf(0))) && ml_exp(ml_make_inf(0)) > 0.0,
        "exp(+inf) is +inf");

    ASSERT_TRUE(&ctx,
        ml_exp(-ml_make_inf(0)) == 0.0,
        "exp(-inf) is 0");

    return ml_test_summary(&ctx);
}
"""


PUNCHLIST_LOG = r"""
<!-- MATHLIB_CLOSURE_P2_P0_3_LOG -->
## Red Team P2 P0-3

- Fixed `ml_exp()` premature overflow threshold.
- Fixed `ml_exp()` premature underflow threshold.
- Added regression suite `tests/test_edge_redteam_p0_3.c`.
"""


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


def insert_macros(v11s: Path, force: bool) -> None:
    path = v11s / "src" / "exp_log.c"

    if not path.is_file():
        fail(f"Missing expected file: {path}")

    original = path.read_text(encoding="utf-8")
    text = normalize(original)

    marker = "MATHLIB_CLOSURE_P2_P0_3_EXP_LIMITS"

    if marker in text:
        print(f"[skip] {path}: P0-3 macro marker already present")
        return

    pattern = re.compile(r'(?m)^(#include "internal/pow_util.h".*)$')

    patched, count = pattern.subn(
        lambda m: m.group(1) + "\n" + MACRO_BLOCK,
        text,
        count=1,
    )

    if count != 1:
        fail(
            f"{path}: could not find internal/pow_util.h include anchor. "
            "Source may have drifted."
        )

    write_text(path, patched)


def patch_exp(v11s: Path, force: bool) -> None:
    path = v11s / "src" / "exp_log.c"

    original = path.read_text(encoding="utf-8")
    text = normalize(original)

    marker = "MATHLIB_CLOSURE_P2_P0_3_EXP_LIMITS"

    if marker in text and not force:
        print(f"[skip] {path}: ml_exp P0-3 marker already present")
        return

    pattern = re.compile(
        r"(?ms)^[ \t]*ML_API[ \t]+double[ \t]+ml_exp\(double[ \t]+x\)[ \t]*\{.*?"
        r"(?=^[ \t]*ML_API[ \t]+double[ \t]+ml_log\(|\Z)"
    )

    patched, count = pattern.subn(
        lambda m: NEW_EXP + "\n",
        text,
        count=1,
    )

    if count != 1:
        fail(
            f"{path}: expected exactly one ml_exp() match, got {count}. "
            "Source may have drifted."
        )

    write_text(path, patched)


def write_test(v11s: Path, force: bool) -> None:
    test_path = v11s / "tests" / "test_edge_redteam_p0_3.c"

    if test_path.exists() and not force:
        try:
            old = test_path.read_text(encoding="utf-8")
            if "MATHLIB_CLOSURE_P2_P0_3_TEST" in old:
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

    if "<!-- MATHLIB_CLOSURE_P2_P0_3_LOG -->" in text:
        print(f"[skip] {punch_path}: P0-3 log already present")
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
    print("  Red Team P2 P0-3: ml_exp threshold fix")
    print("=========================================================")
    print(f"Root: {root}")
    print(f"v11S: {v11s}")
    print(f"Force: {force}")
    print("---------------------------------------------------------")

    insert_macros(v11s, force)
    patch_exp(v11s, force)
    write_test(v11s, force)
    append_log(v11s)
    archive_self(v11s, force)

    print("---------------------------------------------------------")
    print("P0-3 changes applied.")
    print("")
    print("Next verification step:")
    print("")
    print("    cd v11S")
    print("    MATHLIB_EDGE_SANITIZERS=1 bash tests/run_edge_tests.sh")
    print("")
    print("After that passes, the next block is:")
    print("")
    print("    P0-4: fix ml_sinh() / ml_cosh() overflow thresholds")
    print("=========================================================")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
