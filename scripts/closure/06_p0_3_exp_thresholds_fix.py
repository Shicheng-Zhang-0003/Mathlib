#!/usr/bin/env python3
"""
06_p0_3_exp_thresholds_fix.py

Run from the folder that CONTAINS the v11S working folder.

This is the corrected P0-3 script.

It fixes the earlier marker collision by using separate markers for:

    - limit macros
    - ml_exp() function replacement
    - regression test
    - closure log

Targets:

    v11S/src/exp_log.c
    v11S/tests/test_edge_redteam_p0_3.c
    v11S/docs/CLOSURE_PUNCHLIST.md

Usage:

    python3 06_p0_3_exp_thresholds_fix.py

Force overwrite:

    python3 06_p0_3_exp_thresholds_fix.py --force
"""

from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path


MACRO_MARKER = "MATHLIB_CLOSURE_P2_P0_3_EXP_LIMITS"
FUNC_MARKER = "MATHLIB_CLOSURE_P2_P0_3_EXP_FUNC"
TEST_MARKER = "MATHLIB_CLOSURE_P2_P0_3_TEST"
LOG_MARKER = "MATHLIB_CLOSURE_P2_P0_3_FIX_LOG"


MACRO_BLOCK = r"""/* MATHLIB_CLOSURE_P2_P0_3_EXP_LIMITS */
#ifndef ML_LOG_DBL_MAX
#define ML_LOG_DBL_MAX 709.782712893384
#endif

#ifndef ML_LOG_UNDERFLOW
#define ML_LOG_UNDERFLOW (-745.133219101941)
#endif
"""


NEW_EXP = r"""ML_API double ml_exp(double x) {
    /* MATHLIB_CLOSURE_P2_P0_3_EXP_FUNC */
    if (ml_isnan(x)) return x;
    if (ml_isinf(x)) return (x > 0.0) ? ml_make_inf(0) : 0.0;
    if (x == 0.0) return 1.0;

    /*
     * P0-3: use double-limit-aware thresholds.
     *
     * Old code used:
     *
     *   if (x > 709.78) return inf;
     *   if (x < -745.13) return 0;
     *
     * Those thresholds were too conservative.
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
     */
    ASSERT_TRUE(&ctx,
        ml_isfinite(ml_exp(709.781)),
        "exp(709.781) must remain finite");

    ASSERT_TRUE(&ctx,
        ml_isfinite(ml_exp(709.782)),
        "exp(709.782) must remain finite");

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
<!-- MATHLIB_CLOSURE_P2_P0_3_FIX_LOG -->
## Red Team P2 P0-3 Fix

- Corrected earlier P0-3 marker collision.
- `ml_exp()` now uses `ML_LOG_DBL_MAX` and `ML_LOG_UNDERFLOW`.
- Old coarse thresholds `709.78` and `-745.13` are replaced.
- Regression suite `tests/test_edge_redteam_p0_3.c` remains active.
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


def locate_v11s() -> tuple[Path, Path]:
    root = Path.cwd()

    candidate = root / "v11S"
    if candidate.is_dir():
        return root, candidate

    # Convenience: allow running from inside v11S itself.
    if (root / "src" / "exp_log.c").is_file() and (root / "include" / "mathlib").is_dir():
        print("[note] Running from inside v11S; treating current directory as v11S.")
        return root.parent, root

    fail(
        "Run this script from the folder that CONTAINS the v11S directory, "
        "or from inside v11S itself."
    )


def ensure_macros(v11s: Path, force: bool) -> None:
    path = v11s / "src" / "exp_log.c"

    if not path.is_file():
        fail(f"Missing expected file: {path}")

    text = normalize(path.read_text(encoding="utf-8"))

    if MACRO_MARKER in text:
        print(f"[skip] {path}: exp limit macros already present")
        return

    anchors = [
        r'(?m)^(#include "internal/pow_util.h".*)$',
        r'(?m)^(#include "internal/hypot.h".*)$',
        r'(?m)^(#include "ml_exp_log.h".*)$',
    ]

    for anchor in anchors:
        pattern = re.compile(anchor)
        patched, count = pattern.subn(
            lambda m: m.group(1) + "\n" + MACRO_BLOCK,
            text,
            count=1,
        )

        if count == 1:
            write_text(path, patched)
            return

    fail(f"{path}: could not find a safe include anchor for exp limit macros.")


def patch_exp(v11s: Path, force: bool) -> None:
    path = v11s / "src" / "exp_log.c"

    text = normalize(path.read_text(encoding="utf-8"))

    if FUNC_MARKER in text and not force:
        print(f"[skip] {path}: ml_exp P0-3 function marker already present")
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

    if "if (x > ML_LOG_DBL_MAX)" not in patched:
        fail(f"{path}: patch verification failed: missing ML_LOG_DBL_MAX guard.")

    if "if (x < ML_LOG_UNDERFLOW)" not in patched:
        fail(f"{path}: patch verification failed: missing ML_LOG_UNDERFLOW guard.")

    write_text(path, patched)


def write_test(v11s: Path, force: bool) -> None:
    test_path = v11s / "tests" / "test_edge_redteam_p0_3.c"

    if test_path.exists() and not force:
        try:
            old = normalize(test_path.read_text(encoding="utf-8"))
            if TEST_MARKER in old:
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

    if LOG_MARKER in text:
        print(f"[skip] {punch_path}: P0-3 fix log already present")
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

    root, v11s = locate_v11s()

    print("=========================================================")
    print("  Red Team P2 P0-3 Fix: ml_exp threshold correction")
    print("=========================================================")
    print(f"Root: {root}")
    print(f"v11S: {v11s}")
    print(f"Force: {force}")
    print("---------------------------------------------------------")

    ensure_macros(v11s, force)
    patch_exp(v11s, force)
    write_test(v11s, force)
    append_log(v11s)
    archive_self(v11s, force)

    print("---------------------------------------------------------")
    print("P0-3 fix applied.")
    print("")
    print("Quick verification:")
    print("")
    print("    grep -n \"MATHLIB_CLOSURE_P2_P0_3_EXP_FUNC\" v11S/src/exp_log.c")
    print("    grep -n \"x > ML_LOG_DBL_MAX\" v11S/src/exp_log.c")
    print("    grep -n \"x < ML_LOG_UNDERFLOW\" v11S/src/exp_log.c")
    print("")
    print("Then rerun edge tests:")
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
