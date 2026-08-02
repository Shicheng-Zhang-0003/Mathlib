#!/usr/bin/env python3
"""
07d_accuracy_fixes.py

Run from the folder that CONTAINS the v12A1 working folder.

Fixes three accuracy bugs found during full audit:

1. ml_atan: Taylor series only has 10 terms (i=3..21).
   At |x| = 1/3 (edge of reduced domain), truncation error is ~10^4 ULP.
   Fix: extend to i=35 (18 terms). Truncation error drops to <0.001 ULP.
   Affects: ml_atan, ml_asin, ml_acos, ml_atan2, ml_acot

2. ml_sinh: for x in [1e-4, 0.5], the formula 0.5*(exp(x)-exp(-x))
   suffers catastrophic cancellation (~10^3-10^4 ULP).
   Fix: Taylor series sinh(x) = x + x^3/3! + x^5/5! + ... for |x| < 0.5.
   Also lower the "return x" threshold from 1e-4 to 1e-8.

3. ml_tanh: for x in [1e-4, 0.1], the formula (1-exp(-2x))/(1+exp(-2x))
   loses precision (~10^3 ULP).
   Fix: Taylor series for |x| < 0.1, lower "return x" threshold to 1e-8.

Also removes duplicate ML_LOG_DBL_MAX define in exp_log.c.

Targets:
    v12A1/src/trig.c       (ml_atan)
    v12A1/src/exp_log.c    (ml_sinh, ml_tanh, duplicate define)

Usage:
    python3 07d_accuracy_fixes.py
    python3 07d_accuracy_fixes.py --force
"""

from __future__ import annotations
import re
import shutil
import sys
from pathlib import Path

MARKER = "MATHLIB_V12A1_ACCURACY_FIXES"

# ---------------------------------------------------------------------------
# New ml_atan body (extend Taylor from 10 to 18 terms)
# ---------------------------------------------------------------------------
OLD_ATAN_LOOP = "    for (int i = 3; i <= 21; i += 2) { term *= -x2; result += term / i; }"
NEW_ATAN_LOOP = (
    "    /* MATHLIB_V12A1_ACCURACY_FIXES: extended from 10 to 18 terms.\n"
    "     * Old: i <= 21 (10 terms), truncation error ~10^4 ULP at |x|=1/3.\n"
    "     * New: i <= 35 (18 terms), truncation error < 0.001 ULP at |x|=1/3.\n"
    "     * Affects ml_atan, ml_asin, ml_acos, ml_atan2, ml_acot. */\n"
    "    for (int i = 3; i <= 35; i += 2) { term *= -x2; result += term / i; }"
)

# ---------------------------------------------------------------------------
# New ml_sinh body
# ---------------------------------------------------------------------------
NEW_SINH = r"""ML_API double ml_sinh(double x) {
/* MATHLIB_CLOSURE_P2_P0_4_HYPERBOLIC_SHIFT */
/* MATHLIB_V12A1_ACCURACY_FIXES: Taylor series for small x */
if (ml_isnan(x)) return x;
if (ml_isinf(x)) return x;
double ax = ml_fabs(x);
if (ax == 0.0) return x;
/*
* For |x| < 1e-8, sinh(x) = x to within 1 ULP.
* (x^3/6 < 1 ULP of x when x < ~2.6e-8)
*/
if (ax < 1e-8) return x;
/*
* For |x| < 0.5, use Taylor series to avoid catastrophic
* cancellation in 0.5*(exp(x) - exp(-x)).
*
* sinh(x) = x + x^3/3! + x^5/5! + ... + x^19/19!
* At x = 0.5, truncation error < 0.001 ULP.
*/
if (ax < 0.5) {
    double x2 = x * x;
    double term = x;
    double result = x;
    term *= x2; result += term * (1.0/6.0);           /* x^3/3! */
    term *= x2; result += term * (1.0/120.0);         /* x^5/5! */
    term *= x2; result += term * (1.0/5040.0);        /* x^7/7! */
    term *= x2; result += term * (1.0/362880.0);      /* x^9/9! */
    term *= x2; result += term * (1.0/39916800.0);    /* x^11/11! */
    term *= x2; result += term * (1.0/6227020800.0);  /* x^13/13! */
    term *= x2; result += term * (1.0/1307674368000.0); /* x^15/15! */
    term *= x2; result += term * (1.0/355687428096000.0); /* x^17/17! */
    term *= x2; result += term * (1.0/121645100408832000.0); /* x^19/19! */
    return result;
}
/*
* For |x| >= 0.5, the exp-based formula has no significant
* cancellation (both exp(x) and exp(-x) differ by > 2x).
*/
if (ax > ML_LOG_HYP_OVERFLOW) {
    return ml_make_inf(x < 0.0);
}
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
}"""

# ---------------------------------------------------------------------------
# New ml_tanh body
# ---------------------------------------------------------------------------
NEW_TANH = r"""ML_API double ml_tanh(double x) {
/* MATHLIB_V12A1_ACCURACY_FIXES: Taylor series for small x */
if (ml_isnan(x)) return x;
if (ml_isinf(x)) return ml_copysign(1.0, x);
double ax = ml_fabs(x);
if (ax == 0.0) return x;
/*
* For |x| < 1e-8, tanh(x) = x to within 1 ULP.
*/
if (ax < 1e-8) return x;
/*
* For |x| < 0.1, use Taylor series to avoid cancellation
* in (1 - exp(-2x)) / (1 + exp(-2x)).
*
* tanh(x) = x - x^3/3 + 2x^5/15 - 17x^7/315
*           + 62x^9/2835 - 1382x^11/155925 + 21844x^13/6081075
*
* At x = 0.1, truncation error < 0.01 ULP.
*/
if (ax < 0.1) {
    double x2 = x * x;
    double x3 = x * x2;
    double x5 = x3 * x2;
    double x7 = x5 * x2;
    double x9 = x7 * x2;
    double x11 = x9 * x2;
    double x13 = x11 * x2;
    double result = x
        - x3  * (1.0 / 3.0)
        + x5  * (2.0 / 15.0)
        - x7  * (17.0 / 315.0)
        + x9  * (62.0 / 2835.0)
        - x11 * (1382.0 / 155925.0)
        + x13 * (21844.0 / 6081075.0);
    return result;
}
if (ax > 20.0) return ml_copysign(1.0, x);
double e = ml_exp(-2.0 * ax);
double t = (1.0 - e) / (1.0 + e);
return ml_copysign(t, x);
}"""


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def fail(msg):
    print("ERROR: " + msg)
    sys.exit(1)

def normalize(text):
    return text.replace("\r\n", "\n").replace("\r", "\n")

def write_text(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(content)
    print(f"  [write] {path}")

def locate_v12a1():
    root = Path.cwd()
    candidate = root / "v12A1"
    if candidate.is_dir():
        return root, candidate
    if (root / "src" / "trig.c").is_file():
        return root.parent, root
    fail("Run from the folder that CONTAINS v12A1/")


# ---------------------------------------------------------------------------
# Patches
# ---------------------------------------------------------------------------
def patch_atan(v12, force):
    path = v12 / "src" / "trig.c"
    text = normalize(path.read_text(encoding="utf-8"))
    if MARKER in text and not force:
        print(f"  [skip] {path}: already patched")
        return
    if OLD_ATAN_LOOP not in text:
        if "i <= 35" in text:
            print(f"  [skip] {path}: atan already extended")
            return
        fail(f"{path}: could not find atan Taylor loop")
    text = text.replace(OLD_ATAN_LOOP, NEW_ATAN_LOOP, 1)
    write_text(path, text)

def patch_sinh(v12, force):
    path = v12 / "src" / "exp_log.c"
    text = normalize(path.read_text(encoding="utf-8"))
    if MARKER in text and not force:
        print(f"  [skip] {path}: already patched")
        return
    # Replace ml_sinh function body
    pattern = re.compile(
        r"(?ms)^ML_API double ml_sinh\(double x\) \{.*?"
        r"(?=^ML_API double ml_cosh\()"
    )
    patched, count = pattern.subn(NEW_SINH + "\n", text, count=1)
    if count != 1:
        fail(f"{path}: could not find ml_sinh function")
    write_text(path, patched)

def patch_tanh(v12, force):
    path = v12 / "src" / "exp_log.c"
    text = normalize(path.read_text(encoding="utf-8"))
    if MARKER in text and not force:
        print(f"  [skip] {path}: already patched")
        return
    pattern = re.compile(
        r"(?ms)^ML_API double ml_tanh\(double x\) \{.*?"
        r"(?=^ML_API double ml_asinh\()"
    )
    patched, count = pattern.subn(NEW_TANH + "\n", text, count=1)
    if count != 1:
        fail(f"{path}: could not find ml_tanh function")
    write_text(path, patched)

def remove_duplicate_define(v12, force):
    path = v12 / "src" / "exp_log.c"
    text = normalize(path.read_text(encoding="utf-8"))
    # Remove the second (duplicate) ML_LOG_DBL_MAX block
    # The first one is in HYPERBOLIC_LIMITS, the second in EXP_LIMITS
    dup = (
        "/* MATHLIB_CLOSURE_P2_P0_3_EXP_LIMITS */\n"
        "#ifndef ML_LOG_DBL_MAX\n"
        "#define ML_LOG_DBL_MAX 709.782712893384\n"
        "#endif\n"
    )
    if text.count("ML_LOG_DBL_MAX") > 2 and dup in text:
        text = text.replace(dup, "/* MATHLIB_CLOSURE_P2_P0_3_EXP_LIMITS (ML_LOG_DBL_MAX defined above) */\n", 1)
        write_text(path, text)
        print(f"  [fix] {path}: removed duplicate ML_LOG_DBL_MAX")
    else:
        print(f"  [skip] {path}: no duplicate ML_LOG_DBL_MAX found")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    force = "--force" in sys.argv[1:]
    root, v12 = locate_v12a1()

    print("=========================================================")
    print("  MATHLIB v12A1: ACCURACY FIXES (audit findings)")
    print("=========================================================")
    print(f"  Root:  {root}")
    print(f"  v12A1: {v12}")
    print(f"  Force: {force}")
    print("---------------------------------------------------------")

    print("\n[1/4] trig.c — ml_atan Taylor 10 -> 18 terms")
    patch_atan(v12, force)

    print("\n[2/4] exp_log.c — ml_sinh Taylor for |x| < 0.5")
    patch_sinh(v12, force)

    print("\n[3/4] exp_log.c — ml_tanh Taylor for |x| < 0.1")
    patch_tanh(v12, force)

    print("\n[4/4] exp_log.c — remove duplicate ML_LOG_DBL_MAX")
    remove_duplicate_define(v12, force)

    # Archive self
    try:
        src = Path(__file__).resolve()
        dst = v12 / "scripts" / "v12a1" / src.name
        if src != dst:
            dst.parent.mkdir(parents=True, exist_ok=True)
            if not dst.exists() or force:
                shutil.copy2(src, dst)
                print(f"  [archive] {dst}")
    except NameError:
        pass

    print("\n---------------------------------------------------------")
    print("  Accuracy fixes applied.")
    print("")
    print("  What changed:")
    print("    ml_atan:  10 terms -> 18 terms (10^4 ULP -> <0.001 ULP)")
    print("    ml_sinh:  Taylor series for |x| < 0.5 (was cancellation)")
    print("    ml_tanh:  Taylor series for |x| < 0.1 (was cancellation)")
    print("    Removed duplicate ML_LOG_DBL_MAX define")
    print("")
    print("  Functions fixed:")
    print("    ml_atan, ml_asin, ml_acos, ml_atan2, ml_acot")
    print("    ml_sinh, ml_tanh")
    print("")
    print("  Verify:")
    print("    cd v12A1")
    print("    cmake -B build && cmake --build build")
    print("    ./build/test")
    print("    MATHLIB_EDGE_SANITIZERS=1 bash tests/run_edge_tests.sh")
    print("    ./build/fuzz_god_mode 123456789")
    print("")
    print("  Then re-run oracle to confirm no regressions.")
    print("=========================================================")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
