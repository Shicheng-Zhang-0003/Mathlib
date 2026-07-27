#!/usr/bin/env python3
"""
06_p0_4_hyperbolic_thresholds.py

Run from the folder that CONTAINS the v11S working folder.

Expected layout:

    parent_folder/
        06_p0_4_hyperbolic_thresholds.py
        v11S/
            src/
            tests/
            docs/
            ...

This script implements Red Team P2 item P0-4:

    Fix ml_sinh() and ml_cosh() premature overflow thresholds.

It modifies:

    v11S/src/exp_log.c

It adds:

    v11S/tests/test_edge_redteam_p0_4.c

It appends a closure log entry to:

    v11S/docs/CLOSURE_PUNCHLIST.md

Usage:

    python3 06_p0_4_hyperbolic_thresholds.py

Force overwrite:

    python3 06_p0_4_hyperbolic_thresholds.py --force
"""

from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path


MACRO_BLOCK = r"""/* MATHLIB_CLOSURE_P2_P0_4_HYPERBOLIC_LIMITS */
#ifndef ML_LOG_DBL_MAX
#define ML_LOG_DBL_MAX 709.782712893384
#endif

#ifndef ML_LOG_HYP_OVERFLOW
#define ML_LOG_HYP_OVERFLOW (ML_LOG_DBL_MAX + ML_LN2)
#endif
"""


NEW_SINH = r"""ML_API double ml_sinh(double x) {
    /* MATHLIB_CLOSURE_P2_P0_4_HYPERBOLIC_SHIFT */
    if (ml_isnan(x)) return x;
    if (ml_isinf(x)) return x;

    double ax = ml_fabs(x);

    if (ax < 1e-4) {
        return x;
    }

    /*
     * sinh(x) is approximately:
     *
     *   0.5 * exp(x)
     *
     * for large positive x.
     *
     * Therefore overflow happens near:
     *
     *   log(DBL_MAX) + log(2)
     *
     * not merely log(DBL_MAX).
     */
    if (ax > ML_LOG_HYP_OVERFLOW) {
        return ml_make_inf(x < 0.0);
    }

    /*
     * Near overflow, avoid computing exp(ax) directly.
     *
     * Since:
     *
     *   0.5 * exp(ax) = exp(ax - ln2)
     *
     * we can stay finite longer and avoid premature overflow.
     */
    if (ax > 700.0) {
        double ep_half = ml_exp(ax - ML_LN2);
        double em_half = ml_exp(-ax - ML_LN2);
        double r = ep_half - em_half;
        return (x < 0.0) ? -r : r;
    }

    double ep = ml_exp(ax);
    double em = ml_exp(-ax);
    double r = 0.5 * (ep - em);
    return (x < 0.0) ? -r : r;
}
"""


NEW_COSH = r"""ML_API double ml_cosh(double x) {
    /* MATHLIB_CLOSURE_P2_P0_4_HYPERBOLIC_SHIFT */
    if (ml_isnan(x)) return x;
    if (ml_isinf(x)) return ml_make_inf(0);

    double ax = ml_fabs(x);

    /*
     * cosh(x) is approximately:
     *
     *   0.5 * exp(x)
     *
     * for large |x|.
     *
     * Overflow happens near:
     *
     *   log(DBL_MAX) + log(2)
     */
    if (ax > ML_LOG_HYP_OVERFLOW) {
        return ml_make_inf(0);
    }

    /*
     * Near overflow, use the shifted form:
     *
     *   0.5 * exp(ax) = exp(ax - ln2)
     */
    if (ax > 700.0) {
        double ep_half = ml_exp(ax - ML_LN2);
        double em_half = ml_exp(-ax - ML_LN2);
        return ep_half + em_half;
    }

    double ep = ml_exp(ax);
    double em = ml_exp(-ax);
    return 0.5 * (ep + em);
}
"""


TEST_C = r"""/* MATHLIB_CLOSURE_P2_P0_4_TEST */
/* Red Team P2 P0-4: hyperbolic overflow threshold regression tests */
#include "test_harness.h"
#include "ml_core.h"
#include "ml_exp_log.h"

int main(void) {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Edge Red Team P0-4");

    /*
     * Basic special-case behavior.
     */
    ASSERT_TRUE(&ctx,
        ml_isnan(ml_sinh(ml_make_nan())),
        "sinh(NaN) is NaN");

    ASSERT_TRUE(&ctx,
        ml_isnan(ml_cosh(ml_make_nan())),
        "cosh(NaN) is NaN");

    ASSERT_TRUE(&ctx,
        ml_isinf(ml_sinh(ml_make_inf(0))) && ml_sinh(ml_make_inf(0)) > 0.0,
        "sinh(+inf) is +inf");

    ASSERT_TRUE(&ctx,
        ml_isinf(ml_sinh(-ml_make_inf(0))) && ml_sinh(-ml_make_inf(0)) < 0.0,
        "sinh(-inf) is -inf");

    ASSERT_TRUE(&ctx,
        ml_isinf(ml_cosh(ml_make_inf(0))) && ml_cosh(ml_make_inf(0)) > 0.0,
        "cosh(+inf) is +inf");

    ASSERT_TRUE(&ctx,
        ml_isinf(ml_cosh(-ml_make_inf(0))) && ml_cosh(-ml_make_inf(0)) > 0.0,
        "cosh(-inf) is +inf");

    /*
     * The old threshold was approximately log(DBL_MAX):
     *
     *   709.782712893384
     *
     * But sinh/cosh overflow near:
     *
     *   log(DBL_MAX) + log(2)
     *   ~= 710.475860073944
     *
     * Therefore values around 710.0 must remain finite.
     */
    double s710 = ml_sinh(710.0);
    ASSERT_TRUE(&ctx,
        ml_isfinite(s710) && s710 > 0.0,
        "sinh(710) must be finite");

    double c710 = ml_cosh(710.0);
    ASSERT_TRUE(&ctx,
        ml_isfinite(c710) && c710 > 0.0,
        "cosh(710) must be finite");

    double s7104 = ml_sinh(710.4);
    ASSERT_TRUE(&ctx,
        ml_isfinite(s7104) && s7104 > 0.0,
        "sinh(710.4) must be finite");

    double c7104 = ml_cosh(710.4);
    ASSERT_TRUE(&ctx,
        ml_isfinite(c7104) && c7104 > 0.0,
        "cosh(710.4) must be finite");

    /*
     * Beyond the corrected hyperbolic overflow boundary, result must be inf.
     */
    ASSERT_TRUE(&ctx,
        ml_isinf(ml_sinh(710.5)) && ml_sinh(710.5) > 0.0,
        "sinh(710.5) must overflow to +inf");

    ASSERT_TRUE(&ctx,
        ml_isinf(ml_cosh(710.5)) && ml_cosh(710.5) > 0.0,
        "cosh(710.5) must overflow to +inf");

    ASSERT_TRUE(&ctx,
        ml_isinf(ml_sinh(-710.5)) && ml_sinh(-710.5) < 0.0,
        "sinh(-710.5) must overflow to -inf");

    /*
     * Existing large-argument behavior must remain intact.
     */
    ASSERT_TRUE(&ctx,
        ml_isinf(ml_sinh(1000.0)) && ml_sinh(1000.0) > 0.0,
        "sinh(1000) is +inf");

    ASSERT_TRUE(&ctx,
        ml_isinf(ml_sinh(-1000.0)) && ml_sinh(-1000.0) < 0.0,
        "sinh(-1000) is -inf");

    ASSERT_TRUE(&ctx,
        ml_isinf(ml_cosh(1000.0)) && ml_cosh(1000.0) > 0.0,
        "cosh(1000) is +inf");

    ASSERT_TRUE(&ctx,
        ml_isinf(ml_cosh(-1000.0)) && ml_cosh(-1000.0) > 0.0,
        "cosh(-1000) is +inf");

    return ml_test_summary(&ctx);
}
"""


PUNCHLIST_LOG = r"""
<!-- MATHLIB_CLOSURE_P2_P0_4_LOG -->
## Red Team P2 P0-4

- Fixed premature `ml_sinh()` overflow threshold.
- Fixed premature `ml_cosh()` overflow threshold.
- Hyperbolic overflow boundary now uses `log(DBL_MAX) + log(2)`.
- Near-overflow evaluation now uses shifted exponential form `exp(x - ln2)`.
- Added regression suite `tests/test_edge_redteam_p0_4.c`.
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
    if (root / "src" / "exp_log.c").is_file() and (root / "include" / "mathlib").is_dir():
        print("[note] Running from inside v11S; treating current directory as v11S.")
        return root.parent, root

    fail(
        "Run this script from the folder that CONTAINS the v11S directory, "
        "or from inside v11S itself."
    )


def ensure_macros(v11s: Path) -> None:
    path = v11s / "src" / "exp_log.c"

    if not path.is_file():
        fail(f"Missing expected file: {path}")

    text = normalize(path.read_text(encoding="utf-8"))

    marker = "MATHLIB_CLOSURE_P2_P0_4_HYPERBOLIC_LIMITS"

    if marker in text:
        print(f"[skip] {path}: hyperbolic limit macros already present")
        return

    anchor = re.compile(r'(?m)^(#include "internal/pow_util.h".*)$')

    patched, count = anchor.subn(
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


def patch_hyperbolics(v11s: Path, force: bool) -> None:
    path = v11s / "src" / "exp_log.c"

    text = normalize(path.read_text(encoding="utf-8"))

    marker = "MATHLIB_CLOSURE_P2_P0_4_HYPERBOLIC_SHIFT"

    if marker in text and not force:
        print(f"[skip] {path}: ml_sinh/ml_cosh P0-4 marker already present")
        return

    sinh_pattern = re.compile(
        r"(?ms)^[ \t]*ML_API[ \t]+double[ \t]+ml_sinh\(double[ \t]+x\)[ \t]*\{.*?"
        r"(?=^[ \t]*ML_API[ \t]+double[ \t]+ml_cosh\(|\Z)"
    )

    cosh_pattern = re.compile(
        r"(?ms)^[ \t]*ML_API[ \t]+double[ \t]+ml_cosh\(double[ \t]+x\)[ \t]*\{.*?"
        r"(?=^[ \t]*ML_API[ \t]+double[ \t]+ml_tanh\(|\Z)"
    )

    patched, count = sinh_pattern.subn(
        lambda m: NEW_SINH + "\n",
        text,
        count=1,
    )

    if count != 1:
        fail(
            f"{path}: expected exactly one ml_sinh() match, got {count}. "
            "Source may have drifted."
        )

    patched, count = cosh_pattern.subn(
        lambda m: NEW_COSH + "\n",
        patched,
        count=1,
    )

    if count != 1:
        fail(
            f"{path}: expected exactly one ml_cosh() match, got {count}. "
            "Source may have drifted."
        )

    write_text(path, patched)


def write_test(v11s: Path, force: bool) -> None:
    path = v11s / "tests" / "test_edge_redteam_p0_4.c"

    if path.exists() and not force:
        try:
            old = normalize(path.read_text(encoding="utf-8"))
            if "MATHLIB_CLOSURE_P2_P0_4_TEST" in old:
                print(f"[skip] {path}: already present")
                return
        except UnicodeDecodeError:
            pass

    write_text(path, TEST_C)


def append_log(v11s: Path) -> None:
    path = v11s / "docs" / "CLOSURE_PUNCHLIST.md"

    if not path.exists():
        write_text(
            path,
            "# MathLib v11S Closure Punchlist\n" + PUNCHLIST_LOG,
        )
        return

    text = normalize(path.read_text(encoding="utf-8"))

    if "<!-- MATHLIB_CLOSURE_P2_P0_4_LOG -->" in text:
        print(f"[skip] {path}: P0-4 log already present")
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
    print("  Red Team P2 P0-4: hyperbolic overflow threshold fix")
    print("=========================================================")
    print(f"Root: {root}")
    print(f"v11S: {v11s}")
    print(f"Force: {force}")
    print("---------------------------------------------------------")

    ensure_macros(v11s)
    patch_hyperbolics(v11s, force)
    write_test(v11s, force)
    append_log(v11s)
    archive_self(v11s, force)

    print("---------------------------------------------------------")
    print("P0-4 changes applied.")
    print("")
    print("Quick verification:")
    print("")
    print("    grep -n \"MATHLIB_CLOSURE_P2_P0_4_HYPERBOLIC_SHIFT\" v11S/src/exp_log.c")
    print("")
    print("You should see two matches: one in ml_sinh, one in ml_cosh.")
    print("")
    print("Then run edge tests:")
    print("")
    print("    cd v11S")
    print("    MATHLIB_EDGE_SANITIZERS=1 bash tests/run_edge_tests.sh")
    print("")
    print("If using CMake, rebuild fuzzers/tests:")
    print("")
    print("    cmake --build build")
    print("    ./build/fuzz_god_mode 123456789")
    print("")
    print("If using Make:")
    print("")
    print("    make fuzz_god_mode")
    print("    ./build/fuzz_god_mode 123456789")
    print("")
    print("After P0-4 passes, the next block is:")
    print("")
    print("    P0-5: CORDIC / fuzzer non-finite input hardening")
    print("    (or P1, if P0-5 was already applied earlier)")
    print("=========================================================")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
