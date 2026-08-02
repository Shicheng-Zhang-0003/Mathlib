#!/usr/bin/env python3
"""
08_gamma_lanczos.py

Run from the folder that CONTAINS the v12A1 working folder.

Replaces the degree-8 polynomial gamma sketch with a Lanczos
approximation (g=7, n=9). Adds reflection formula for negative
non-integer arguments. Adds ml_lgamma as a new API.

Targets:
    v12A1/src/integral.c
    v12A1/include/mathlib/ml_integral.h
    v12A1/tests/test_edge_integral.c
    v12A1/tests/test.c

Usage:
    python3 08_gamma_lanczos.py
    python3 08_gamma_lanczos.py --force
"""

from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path

MARKER = "MATHLIB_V12A1_GAMMA_LANCZOS"

# ---------------------------------------------------------------------------
# New C code: replaces everything from ml_gamma_new to end of file
# ---------------------------------------------------------------------------
NEW_GAMMA_SECTION = r"""/* MATHLIB_V12A1_GAMMA_LANCZOS */
/*
 * Lanczos approximation for the gamma function.
 *
 * Uses g=7, n=9 coefficients (Numerical Recipes / Godfrey).
 * Achieves ~15 significant digits across the positive real line.
 *
 * The old v11S implementation used a degree-8 polynomial on [1,2]
 * with iterative reduction, giving only ~1e-2 relative accuracy.
 * This is a 13-order-of-magnitude improvement.
 */

static const double ml_lanczos_coeff[9] = {
     0.99999999999980993,
   676.5203681218851,
  -1259.1392167224028,
   771.32342877765313,
  -176.61502916214059,
    12.507343278686905,
    -0.13857109526572012,
     9.9843695780195716e-6,
     1.5056327351493116e-7
};

#define ML_LANCZOS_G 7.0
#define ML_GAMMA_OVERFLOW 171.6243769563027

/*
 * Internal: log Gamma(x) via Lanczos for x > 0.
 *
 * log Gamma(x) = 0.5*log(2pi) + (z+0.5)*log(t) - t + log(Ag(z))
 * where z = x - 1, t = z + g + 0.5.
 */
static double ml_lgamma_positive(double x) {
    double z = x - 1.0;
    double ag = ml_lanczos_coeff[0];
    for (int i = 1; i < 9; i++) {
        ag += ml_lanczos_coeff[i] / (z + (double)i);
    }
    double t = z + ML_LANCZOS_G + 0.5;
    return 0.91893853320467274178  /* 0.5 * log(2*pi) */
         + (z + 0.5) * ml_log(t)
         - t
         + ml_log(ag);
}

ML_API double ml_lgamma(double x) {
    /* MATHLIB_V12A1_LGAMMA */
    if (ml_isnan(x)) return x;
    if (ml_isinf(x)) return ml_make_inf(0);

    /* Poles at zero and negative integers */
    if (x <= 0.0 && x == ml_round(x)) return ml_make_inf(0);

    if (x > 0.0) {
        return ml_lgamma_positive(x);
    }

    /* Reflection: log|Gamma(x)| = log(pi) - log|sin(pi*x)| - log Gamma(1-x) */
    double sinpx = ml_sin(ML_PI * x);
    if (sinpx == 0.0) return ml_make_inf(0);
    return ml_log(ML_PI / ml_fabs(sinpx)) - ml_lgamma_positive(1.0 - x);
}

ML_API double ml_gamma_new(double x) {
    /* MATHLIB_V12A1_GAMMA_LANCZOS */
    if (ml_isnan(x)) return x;
    if (ml_isinf(x)) return x > 0.0 ? ml_make_inf(0) : ml_make_nan();

    /* Poles at zero and negative integers */
    if (x <= 0.0 && x == ml_round(x)) return ml_make_nan();

    if (x > 0.0) {
        if (x > ML_GAMMA_OVERFLOW) return ml_make_inf(0);
        return ml_exp(ml_lgamma_positive(x));
    }

    /* Reflection: Gamma(x) = pi / (sin(pi*x) * Gamma(1-x)) */
    double sinpx = ml_sin(ML_PI * x);
    if (sinpx == 0.0) return ml_make_nan();
    if (1.0 - x > ML_GAMMA_OVERFLOW) {
        /* Gamma(1-x) overflows, so Gamma(x) -> 0 */
        return ml_copysign(0.0, sinpx);
    }
    return ML_PI / (sinpx * ml_exp(ml_lgamma_positive(1.0 - x)));
}
"""

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def fail(message: str) -> None:
    print("ERROR: " + message)
    sys.exit(1)


def normalize(text: str) -> str:
    return text.replace("\r\n", "\n").replace("\r", "\n")


def write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(content)
    print(f"  [write] {path}")


def locate_v12a1() -> tuple[Path, Path]:
    root = Path.cwd()
    candidate = root / "v12A1"
    if candidate.is_dir():
        return root, candidate
    if (root / "src" / "integral.c").is_file():
        print("  [note] Running from inside v12A1.")
        return root.parent, root
    fail("Run from the folder that CONTAINS v12A1/, or from inside v12A1/ itself.")


# ---------------------------------------------------------------------------
# Patches
# ---------------------------------------------------------------------------
def patch_integral(v12: Path, force: bool) -> None:
    path = v12 / "src" / "integral.c"
    if not path.is_file():
        fail(f"Missing expected file: {path}")
    text = normalize(path.read_text(encoding="utf-8"))

    if MARKER in text and not force:
        print(f"  [skip] {path}: already patched")
        return

    # Add ml_trig.h include (needed for reflection formula)
    if '#include "ml_trig.h"' not in text:
        text = text.replace(
            '#include "ml_integral.h"',
            '#include "ml_integral.h"\n#include "ml_trig.h"',
            1,
        )

    # Replace everything from ml_gamma_new to end of file
    pattern = re.compile(
        r"(?ms)^[ \t]*/\*.*?v11S contract:.*?ml_gamma_new.*?\*/\s*"
        r"^[ \t]*ML_API[ \t]+double[ \t]+ml_gamma_new\(double[ \t]+x\)"
        r"[ \t]*\{.*\Z"
    )
    patched, count = pattern.subn(
        lambda m: NEW_GAMMA_SECTION,
        text,
        count=1,
    )
    if count != 1:
        # Fallback: match just the function signature to end of file
        pattern2 = re.compile(
            r"(?ms)^[ \t]*ML_API[ \t]+double[ \t]+ml_gamma_new"
            r"\(double[ \t]+x\)[ \t]*\{.*\Z"
        )
        patched, count = pattern2.subn(
            lambda m: NEW_GAMMA_SECTION,
            text,
            count=1,
        )
        if count != 1:
            fail(f"{path}: could not find ml_gamma_new. Source may have drifted.")

    write_text(path, patched)


def patch_header(v12: Path, force: bool) -> None:
    path = v12 / "include" / "mathlib" / "ml_integral.h"
    if not path.is_file():
        fail(f"Missing expected file: {path}")
    text = normalize(path.read_text(encoding="utf-8"))

    if "ml_lgamma" in text and not force:
        print(f"  [skip] {path}: ml_lgamma already declared")
        return

    text = text.replace(
        "ML_API double ml_gamma_new(double x);",
        "ML_API double ml_gamma_new(double x);\nML_API double ml_lgamma(double x);",
        1,
    )
    write_text(path, text)


def patch_edge_test(v12: Path, force: bool) -> None:
    path = v12 / "tests" / "test_edge_integral.c"
    if not path.is_file():
        fail(f"Missing expected file: {path}")
    text = normalize(path.read_text(encoding="utf-8"))

    if "MATHLIB_V12A1_GAMMA_LANCZOS_TEST" in text and not force:
        print(f"  [skip] {path}: already patched")
        return

    # Tighten existing gamma tolerances
    text = text.replace(
        'ASSERT_NEAR(&ctx, ml_gamma_new(1.0), 1.0, 1e-12, "gamma(1)");',
        'ASSERT_NEAR(&ctx, ml_gamma_new(1.0), 1.0, 1e-14, "gamma(1)");',
    )
    text = text.replace(
        'ASSERT_NEAR(&ctx, ml_gamma_new(2.0), 1.0, 1e-12, "gamma(2)");',
        'ASSERT_NEAR(&ctx, ml_gamma_new(2.0), 1.0, 1e-14, "gamma(2)");',
    )

    # Add new tests before the final return
    new_tests = """
    /* MATHLIB_V12A1_GAMMA_LANCZOS_TEST */
    /* Lanczos accuracy: known exact values */
    ASSERT_NEAR(&ctx, ml_gamma_new(3.0), 2.0, 1e-13, "gamma(3) == 2");
    ASSERT_NEAR(&ctx, ml_gamma_new(4.0), 6.0, 1e-12, "gamma(4) == 6");
    ASSERT_NEAR(&ctx, ml_gamma_new(5.0), 24.0, 1e-11, "gamma(5) == 24");
    ASSERT_NEAR(&ctx, ml_gamma_new(6.0), 120.0, 1e-10, "gamma(6) == 120");
    ASSERT_NEAR(&ctx, ml_gamma_new(0.5), 1.7724538509055159, 1e-13, "gamma(0.5) == sqrt(pi)");

    /* Reflection formula: negative non-integer arguments */
    ASSERT_NEAR(&ctx, ml_gamma_new(-0.5), -3.5449077018110318, 1e-12, "gamma(-0.5)");
    ASSERT_NEAR(&ctx, ml_gamma_new(-1.5), 2.3632718012073548, 1e-12, "gamma(-1.5)");
    ASSERT_TRUE(&ctx, ml_isnan(ml_gamma_new(-1.0)), "gamma(-1) pole is NaN");
    ASSERT_TRUE(&ctx, ml_isnan(ml_gamma_new(-2.0)), "gamma(-2) pole is NaN");

    /* lgamma basic checks */
    ASSERT_NEAR(&ctx, ml_lgamma(1.0), 0.0, 1e-14, "lgamma(1) == 0");
    ASSERT_NEAR(&ctx, ml_lgamma(2.0), 0.0, 1e-14, "lgamma(2) == 0");
    ASSERT_NEAR(&ctx, ml_lgamma(5.0), ml_log(24.0), 1e-12, "lgamma(5) == log(24)");
    ASSERT_TRUE(&ctx, ml_isinf(ml_lgamma(0.0)), "lgamma(0) is inf");
    ASSERT_TRUE(&ctx, ml_isnan(ml_lgamma(ml_make_nan())), "lgamma(NaN) is NaN");

"""
    text = text.replace(
        "    return ml_test_summary(&ctx);",
        new_tests + "    return ml_test_summary(&ctx);",
        1,
    )
    write_text(path, text)


def patch_smoke_test(v12: Path, force: bool) -> None:
    path = v12 / "tests" / "test.c"
    if not path.is_file():
        fail(f"Missing expected file: {path}")
    text = normalize(path.read_text(encoding="utf-8"))

    if "MATHLIB_V12A1_GAMMA_TIGHT" in text and not force:
        print(f"  [skip] {path}: already patched")
        return

    # Tighten gamma tolerances in monolithic smoke test
    text = text.replace(
        'ASSERT_NEAR(&ctx, ml_gamma_new(1.0), 1.0, 1e-3, "gamma(1)");',
        '/* MATHLIB_V12A1_GAMMA_TIGHT */\n'
        '    ASSERT_NEAR(&ctx, ml_gamma_new(1.0), 1.0, 1e-13, "gamma(1)");',
    )
    text = text.replace(
        'ASSERT_NEAR(&ctx, ml_gamma_new(4.0), 6.0, 1e-2, "gamma(4)");',
        'ASSERT_NEAR(&ctx, ml_gamma_new(4.0), 6.0, 1e-11, "gamma(4)");',
    )
    text = text.replace(
        'ASSERT_NEAR(&ctx, ml_gamma_new(5.0), 24.0, 5e-2, "gamma(5)");',
        'ASSERT_NEAR(&ctx, ml_gamma_new(5.0), 24.0, 1e-10, "gamma(5)");',
    )
    write_text(path, text)


def archive_self(v12: Path, force: bool) -> None:
    try:
        source = Path(__file__).resolve()
        dest = v12 / "scripts" / "v12a1" / source.name
        if source == dest:
            return
        if dest.exists() and not force:
            print(f"  [skip] {dest}: already archived")
            return
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, dest)
        print(f"  [archive] {dest}")
    except NameError:
        pass


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    force = "--force" in sys.argv[1:]
    root, v12 = locate_v12a1()

    print("=========================================================")
    print("  MATHLIB v12A1: GAMMA FUNCTION REDESIGN (LANCZOS)")
    print("=========================================================")
    print(f"  Root:  {root}")
    print(f"  v12A1: {v12}")
    print(f"  Force: {force}")
    print("---------------------------------------------------------")

    print("\n[1/4] integral.c — Lanczos gamma + lgamma")
    patch_integral(v12, force)

    print("\n[2/4] ml_integral.h — ml_lgamma declaration")
    patch_header(v12, force)

    print("\n[3/4] test_edge_integral.c — tightened + new tests")
    patch_edge_test(v12, force)

    print("\n[4/4] test.c — smoke test tolerances")
    patch_smoke_test(v12, force)

    archive_self(v12, force)

    print("\n---------------------------------------------------------")
    print("  Gamma redesign applied.")
    print("")
    print("  What changed:")
    print("    - ml_gamma_new: degree-8 poly -> Lanczos g=7 n=9")
    print("    - Reflection formula for negative non-integers")
    print("    - New API: ml_lgamma (log-gamma)")
    print("    - Test tolerances tightened from 1e-2 to 1e-10+")
    print("")
    print("  Verify:")
    print("    cd v12A1")
    print("    cmake -B build && cmake --build build")
    print("    ./build/test")
    print("    MATHLIB_EDGE_SANITIZERS=1 bash tests/run_edge_tests.sh")
    print("")
    print("  Next: 07_payne_hanek.py")
    print("=========================================================")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
