#!/usr/bin/env python3
"""
22_gamma_dd_rewrite.py
Run from the folder that CONTAINS the v12A1 working folder.

GIANT ERROR #5 (second attempt, full rewrite).

Script 21's double-double rewrite failed the oracle:
  - gamma(50/100/171) still 141/316/582 ULP. gamma(171) at 582 ULP
    exceeds even the 257-ULP ceiling of "perfect lgamma rounded to
    double before exp", so the DD low bits never reached ml_exp.
  - Several small-x values REGRESSED (gamma(0.01) 9 -> 18 ULP,
    lgamma(-2.5) 10 -> 140 ULP), which a correct DD chain cannot do.
Conclusion: script 21's DD machinery was buggy. Patching on top is
unsafe, so this script replaces src/integral.c IN FULL.

Design:
  - Full double-double primitive set (two-sum, renorm, DD+DD, DD*DD,
    DD*double, double/double with one Newton correction).
  - ml_log_dd(): true ~106-bit logarithm.
      * z = (m-1)/(m+1): numerator exact (Sterbenz), denominator exact
        via two-sum, quotient DD via one Newton correction.
      * atanh series evaluated with DD Horner.
      * e*ln2: rounding of e*LN2_HI captured by FMA, plus e*LN2_LO.
  - Lanczos Ag sum in DD with DD term divisions.
  - L = 0.5*ln(2pi) - t + w*log(t) + log(Ag) assembled in DD.
  - gamma(x) = ml_exp(L.hi) * (1 + L.lo) via one FMA, so the low part
    actually reaches the exponential (script 21 lost it).
  - Recurrence gamma(x) = gamma(x+1)/x for x < 0.5 (also prevents the
    c1/(z+1) division-by-near-zero catastrophe for tiny x).
  - Exact integer shortcuts: gamma(n) = (n-1)! for n <= 23 (all such
    factorials are exactly representable), lgamma(1) = lgamma(2) = 0,
    lgamma(n) = log((n-1)!) for n <= 23.

Targets:
  v12A1/src/integral.c           (full rewrite)
  v12A1/tests/test_edge_integral.c (full rewrite)

Usage:
  python3 22_gamma_dd_rewrite.py
  python3 22_gamma_dd_rewrite.py --force
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

MARKER = "MATHLIB_V12A1_GAMMA_DD2"
TEST_MARKER = "MATHLIB_V12A1_GAMMA_DD2_TEST"

NEW_INTEGRAL_C = r"""#include "ml_compiler.h"
#include "ml_integral.h"
#include "ml_trig.h"

/*
 * v11S contract:
 * ml_factorial_float(x) means x! = Gamma(x + 1).
 */
ML_API double ml_factorial_float(double x) {
    if (ml_isnan(x)) return x;
    if (x < 0.0) return ml_make_nan();
    if (ml_isinf(x)) return ml_make_inf(0);
    if (x == 0.0) return 1.0;
    return ml_gamma_new(x + 1.0);
}

/*
 * Educational / experimental numerical integrator.
 * Retained for compatibility.
 */
ML_API double ml_integral_traditional(double a, double b, double exponent, double additive, double d) {
    if (ml_isnan(a) || ml_isnan(b) || ml_isnan(exponent) ||
        ml_isnan(additive) || ml_isnan(d)) {
        return ml_make_nan();
    }
    if (ml_isinf(a) || ml_isinf(b) || ml_isinf(exponent) ||
        ml_isinf(additive) || ml_isinf(d)) {
        return ml_make_nan();
    }
    if (d == 0.0) {
        return ml_make_nan();
    }
    if ((d > 0.0 && a >= b) || (d < 0.0 && a <= b)) {
        return 0.0;
    }
    double result = 0.0;
    double x = a;
    const int max_steps = 10000000;
    for (int step = 0; step < max_steps; step++) {
        if ((d > 0.0 && x >= b) || (d < 0.0 && x <= b)) {
            return result;
        }
        double term = ml_pow(x, exponent) + additive;
        if (ml_isnan(term)) {
            return ml_make_nan();
        }
        result += term * d;
        double next_x = x + d;
        if (next_x == x) {
            return ml_make_nan();
        }
        x = next_x;
    }
    return ml_make_nan();
}

/* MATHLIB_V12A1_GAMMA_DD2 */
/*
 * Gamma / log-gamma via Lanczos (g=7, n=9), evaluated end-to-end in
 * double-double arithmetic.
 *
 * Why double-double:
 *   lgamma(x) = 0.5*ln(2pi) + (x-0.5)*ln(x+6.5) - (x+6.5) + ln(Ag(x-1))
 *
 * The product (x-0.5)*ln(x+6.5) amplifies any error in ln by (x-0.5),
 * up to ~170.5. A 53-bit ln(x) carries ~1e-16 relative error, which
 * becomes ~1.7e-14 absolute in the product near x=171 -- hundreds of
 * ULP in gamma(x). Carrying ln and the whole sum in double-double
 * (~106 bits) removes the amplification entirely.
 *
 * The low part must also survive into ml_exp:
 *   gamma = exp(L.hi) * (1 + L.lo)   (one FMA)
 * Rounding L to a single double before exp would throw away up to
 * 0.5 ULP of lgamma, i.e. up to ~257 ULP of gamma at x=171.
 */

#ifndef ML_LN2_HI
#define ML_LN2_HI 6.93147180369123816490e-01
#endif
#ifndef ML_LN2_LO
#define ML_LN2_LO 1.90821492927058500170e-10
#endif

#define ML_LANCZOS_G 7.0
#define ML_GAMMA_OVERFLOW 171.6243769563027
#define ML_HALF_LOG_2PI 0.91893853320467274178

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

/* ------------------------------------------------------------------ */
/* Double-double primitives: a value is represented as hi + lo.        */
/* ------------------------------------------------------------------ */

typedef struct { double hi; double lo; } ml_dd_t;

static inline ml_dd_t ml_dd_from_d(double a) {
    ml_dd_t r;
    r.hi = a;
    r.lo = 0.0;
    return r;
}

/* Exact a + b (Knuth two-sum). */
static inline ml_dd_t ml_dd_two_sum(double a, double b) {
    double s = a + b;
    double v = s - a;
    ml_dd_t r;
    r.hi = s;
    r.lo = (a - (s - v)) + (b - v);
    return r;
}

/* Renormalize hi + lo into canonical double-double form. */
static inline ml_dd_t ml_dd_renorm(double hi, double lo) {
    double s = hi + lo;
    ml_dd_t r;
    r.hi = s;
    r.lo = lo - (s - hi);
    return r;
}

static inline ml_dd_t ml_dd_add(ml_dd_t a, ml_dd_t b) {
    ml_dd_t s = ml_dd_two_sum(a.hi, b.hi);
    return ml_dd_renorm(s.hi, s.lo + a.lo + b.lo);
}

static inline ml_dd_t ml_dd_sub(ml_dd_t a, ml_dd_t b) {
    ml_dd_t nb;
    nb.hi = -b.hi;
    nb.lo = -b.lo;
    return ml_dd_add(a, nb);
}

/* DD * double. Error ~ 2^-106 relative. */
static inline ml_dd_t ml_dd_mul_d(ml_dd_t a, double b) {
    double p = a.hi * b;
    double e = ML_FMA(a.hi, b, -p) + a.lo * b;
    return ml_dd_renorm(p, e);
}

/* DD * DD. Error ~ 2^-106 relative. */
static inline ml_dd_t ml_dd_mul(ml_dd_t a, ml_dd_t b) {
    double p = a.hi * b.hi;
    double e = ML_FMA(a.hi, b.hi, -p) + (a.hi * b.lo + a.lo * b.hi);
    return ml_dd_renorm(p, e);
}

/* double / double -> DD. One Newton correction; ~2^-105 relative. */
static inline ml_dd_t ml_dd_div(double a, double b) {
    double q = a / b;
    double r = ML_FMA(-q, b, a); /* a - q*b; exact (Sterbenz) */
    double q2 = r / b;
    return ml_dd_renorm(q, q2);
}

/* ------------------------------------------------------------------ */
/* ~106-bit natural logarithm for x > 0, finite.                      */
/* ------------------------------------------------------------------ */

static ml_dd_t ml_log_dd(double x) {
    int e;
    double m = ml_frexp_pure(x, &e);
    int adjust = (m < 0.7071067811865475);
    m *= (1.0 + (double)adjust); /* exact scaling into [sqrt(1/2), sqrt(2)] */
    e -= adjust;

    /* z = (m-1)/(m+1).
     * Numerator: exact by Sterbenz (m in [0.5, 2]).
     * Denominator: exact double-double via two-sum.
     * Quotient: one Newton correction gives ~2^-105 relative. */
    double num = m - 1.0;
    ml_dd_t den = ml_dd_two_sum(m, 1.0);
    double q = num / den.hi;
    double r = ML_FMA(-q, den.hi, num);
    r = ML_FMA(-q, den.lo, r);
    double q2 = r / den.hi;
    ml_dd_t z = ml_dd_renorm(q, q2);

    ml_dd_t z2 = ml_dd_mul(z, z);

    /* log(m) = z * P(z^2), P(u) = sum_{k=0..10} (2/(2k+1)) u^k.
     * Truncation error < 2.5e-19 on |z| <= 0.1716. DD Horner keeps
     * evaluation error at ~1e-29. */
    static const double lc[11] = {
        2.0,
        0.6666666666666666,
        0.4,
        0.2857142857142857,
        0.2222222222222222,
        0.18181818181818182,
        0.15384615384615385,
        0.13333333333333333,
        0.11764705882352941,
        0.10526315789473684,
        0.09523809523809523
    };
    ml_dd_t p = ml_dd_from_d(lc[10]);
    for (int i = 9; i >= 0; i--) {
        p = ml_dd_add(ml_dd_mul(p, z2), ml_dd_from_d(lc[i]));
    }
    ml_dd_t lm = ml_dd_mul(z, p);

    /* e*ln2 in DD: e*LN2_HI rounds once; the FMA captures that rounding
     * exactly, and e*LN2_LO restores the truncated tail of ln2. */
    double ed = (double)e;
    double ehi = ed * ML_LN2_HI;
    double elo = ML_FMA(ed, ML_LN2_HI, -ehi) + ed * ML_LN2_LO;
    ml_dd_t eln2 = ml_dd_renorm(ehi, elo);

    return ml_dd_add(lm, eln2);
}

/* ------------------------------------------------------------------ */
/* Exact small factorials and the DD Lanczos core.                    */
/* ------------------------------------------------------------------ */

/* (n-1)! for 1 <= n <= 23. Every intermediate factorial up to 22! is
 * exactly representable in double, so the product is exact. */
static double ml_factorial_exact_small(int n) {
    double f = 1.0;
    for (int k = 2; k <= n - 1; k++) {
        f *= (double)k;
    }
    return f;
}

/* log Gamma(x) for x > 0 as double-double. */
static ml_dd_t ml_lgamma_lanczos_dd(double x) {
    double z = x - 1.0;                    /* exact for x < 2^53 */
    double t = z + ML_LANCZOS_G + 0.5;     /* exact for x < 2^51 */
    double w = z + 0.5;                    /* exact for x < 2^52 */

    /* Ag = c0 + sum c_i/(z+i), accumulated in DD. */
    ml_dd_t ag = ml_dd_from_d(ml_lanczos_coeff[0]);
    for (int i = 1; i < 9; i++) {
        ml_dd_t term = ml_dd_div(ml_lanczos_coeff[i], z + (double)i);
        ag = ml_dd_add(ag, term);
    }
    if (ag.hi <= 0.0) {
        return ml_dd_from_d(ml_make_nan());
    }

    /* log(Ag) = log(ag.hi) + ag.lo/ag.hi + O((ag.lo/ag.hi)^2) */
    ml_dd_t lag = ml_log_dd(ag.hi);
    lag = ml_dd_add(lag, ml_dd_from_d(ag.lo / ag.hi));

    ml_dd_t lt = ml_log_dd(t);
    ml_dd_t wlt = ml_dd_mul_d(lt, w);

    /* L = 0.5*ln(2pi) - t + w*log(t) + log(Ag) */
    ml_dd_t L = ml_dd_two_sum(ML_HALF_LOG_2PI, -t);
    L = ml_dd_add(L, wlt);
    L = ml_dd_add(L, lag);
    return L;
}

/* ------------------------------------------------------------------ */
/* Public APIs                                                         */
/* ------------------------------------------------------------------ */

ML_API double ml_lgamma(double x) {
    if (ml_isnan(x)) return x;
    if (ml_isinf(x)) return ml_make_inf(0);
    /* Poles at zero and negative integers */
    if (x <= 0.0 && x == ml_round(x)) return ml_make_inf(0);
    /* Exact zeros */
    if (x == 1.0 || x == 2.0) return 0.0;

    if (x > 0.0) {
        /* lgamma(n) = log((n-1)!) via the exact factorial */
        if (x == ml_round(x) && x <= 23.0) {
            return ml_log(ml_factorial_exact_small((int)x));
        }
        ml_dd_t L;
        if (x < 0.5) {
            /* lgamma(x) = lgamma(x+1) - log(x) */
            L = ml_lgamma_lanczos_dd(x + 1.0);
            ml_dd_t r = ml_dd_sub(L, ml_dd_from_d(ml_log(x)));
            return r.hi + r.lo;
        }
        L = ml_lgamma_lanczos_dd(x);
        return L.hi + L.lo;
    }

    /* Reflection: log|Gamma(x)| = log(pi/|sin(pi x)|) - lgamma(1-x) */
    double sinpx = ml_sin(ML_PI * x);
    if (sinpx == 0.0) return ml_make_inf(0);
    double logterm = ml_log(ML_PI / ml_fabs(sinpx));
    ml_dd_t L = ml_lgamma_lanczos_dd(1.0 - x);
    ml_dd_t r = ml_dd_sub(ml_dd_from_d(logterm), L);
    return r.hi + r.lo;
}

ML_API double ml_gamma_new(double x) {
    if (ml_isnan(x)) return x;
    if (ml_isinf(x)) return x > 0.0 ? ml_make_inf(0) : ml_make_nan();
    /* Poles at zero and negative integers */
    if (x <= 0.0 && x == ml_round(x)) return ml_make_nan();

    if (x > 0.0) {
        if (x > ML_GAMMA_OVERFLOW) return ml_make_inf(0);
        /* Exact factorial shortcut */
        if (x == ml_round(x) && x <= 23.0) {
            return ml_factorial_exact_small((int)x);
        }
        ml_dd_t L;
        if (x < 0.5) {
            /* Gamma(x) = Gamma(x+1)/x. Also keeps z+1 = x+1 >= 1,
             * avoiding c1/(z+1) -> c1/0 for tiny x. */
            L = ml_lgamma_lanczos_dd(x + 1.0);
            double g = ml_exp(L.hi);
            double G = ML_FMA(g, L.lo, g);
            return G / x;
        }
        L = ml_lgamma_lanczos_dd(x);
        /* gamma = exp(L.hi) * (1 + L.lo): the DD low part reaches the
         * exponential instead of being rounded away. */
        double g = ml_exp(L.hi);
        return ML_FMA(g, L.lo, g);
    }

    /* Reflection: Gamma(x) = pi / (sin(pi x) * Gamma(1-x)) */
    double sinpx = ml_sin(ML_PI * x);
    if (sinpx == 0.0) return ml_make_nan();
    if (1.0 - x > ML_GAMMA_OVERFLOW) {
        /* Gamma(1-x) overflows, so Gamma(x) -> 0 */
        return ml_copysign(0.0, sinpx);
    }
    ml_dd_t L = ml_lgamma_lanczos_dd(1.0 - x);
    double g = ml_exp(L.hi);
    double G = ML_FMA(g, L.lo, g);
    return ML_PI / (sinpx * G);
}
"""

NEW_EDGE_TEST_C = r"""/* v11S CLOSURE IP-20: edge integral tests */
#include "test_harness.h"
#include "ml_integral.h"

int main(void) {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Edge Integral");

    ASSERT_NEAR(&ctx, ml_factorial_float(0.0), 1.0, 1e-15, "factorial_float(0)");
    ASSERT_NEAR(&ctx, ml_factorial_float(5.0), 120.0, 1e-9, "factorial_float(5)");
    ASSERT_TRUE(&ctx, ml_isnan(ml_factorial_float(-1.0)), "factorial_float negative is NaN");

    ASSERT_NEAR(&ctx, ml_gamma_new(1.0), 1.0, 1e-14, "gamma(1)");
    ASSERT_NEAR(&ctx, ml_gamma_new(2.0), 1.0, 1e-14, "gamma(2)");
    ASSERT_TRUE(&ctx, ml_isnan(ml_gamma_new(0.0)), "gamma(0) is NaN");
    double gi = ml_gamma_new(ml_make_inf(0));
    ASSERT_TRUE(&ctx, ml_isinf(gi) && gi > 0.0, "gamma(+inf) is +inf");

    ASSERT_TRUE(&ctx, ml_isnan(ml_integral_traditional(0.0, 1.0, 2.0, 0.0, 0.0)), "integral d=0 is NaN");
    ASSERT_NEAR(&ctx, ml_integral_traditional(0.0, 1.0, 2.0, 0.0, 0.0001), 1.0 / 3.0, 1e-3, "integral x^2");

    /* MATHLIB_V12A1_GAMMA_DD2_TEST */
    /* Exact integer factorials (shortcut path) */
    ASSERT_TRUE(&ctx, ml_gamma_new(3.0) == 2.0, "gamma(3) == 2 exact");
    ASSERT_TRUE(&ctx, ml_gamma_new(4.0) == 6.0, "gamma(4) == 6 exact");
    ASSERT_TRUE(&ctx, ml_gamma_new(5.0) == 24.0, "gamma(5) == 24 exact");
    ASSERT_TRUE(&ctx, ml_gamma_new(6.0) == 120.0, "gamma(6) == 120 exact");
    ASSERT_TRUE(&ctx, ml_gamma_new(10.0) == 362880.0, "gamma(10) == 9! exact");
    ASSERT_TRUE(&ctx, ml_gamma_new(20.0) == 121645100408832000.0, "gamma(20) == 19! exact");
    ASSERT_TRUE(&ctx, ml_gamma_new(23.0) == 1.12400072777760768e21, "gamma(23) == 22! exact");

    /* Tight comparisons against mpmath oracle values (~5.4 ULP tolerance) */
    ASSERT_NEAR(&ctx, ml_gamma_new(0.5), 1.77245385090551610e+00,
                1.77245385090551610e+00 * 1.2e-15, "gamma(0.5) == sqrt(pi) tight");
    ASSERT_NEAR(&ctx, ml_gamma_new(0.1), 9.51350769866873058e+00,
                9.51350769866873058e+00 * 1.2e-15, "gamma(0.1) tight");
    ASSERT_NEAR(&ctx, ml_gamma_new(0.01), 9.94325851191506018e+01,
                9.94325851191506018e+01 * 1.2e-15, "gamma(0.01) tight");
    ASSERT_NEAR(&ctx, ml_gamma_new(0.001), 9.99423772484595474e+02,
                9.99423772484595474e+02 * 1.2e-15, "gamma(0.001) tight");
    ASSERT_NEAR(&ctx, ml_gamma_new(50.0), 6.08281864034267522e+62,
                6.08281864034267522e+62 * 1.2e-15, "gamma(50) tight");
    ASSERT_NEAR(&ctx, ml_gamma_new(100.0), 9.33262154439441533e+155,
                9.33262154439441533e+155 * 1.2e-15, "gamma(100) tight");
    ASSERT_NEAR(&ctx, ml_gamma_new(171.0), 7.25741561530799904e+306,
                7.25741561530799904e+306 * 1.2e-15, "gamma(171) tight");

    /* Reflection formula: negative non-integer arguments */
    ASSERT_NEAR(&ctx, ml_gamma_new(-0.5), -3.54490770181103221e+00,
                3.54490770181103221e+00 * 1.2e-15, "gamma(-0.5) tight");
    ASSERT_NEAR(&ctx, ml_gamma_new(-1.5), 2.36327180120735481e+00,
                2.36327180120735481e+00 * 1.2e-15, "gamma(-1.5) tight");
    ASSERT_NEAR(&ctx, ml_gamma_new(-2.5), -9.45308720482941900e-01,
                9.45308720482941900e-01 * 1.2e-15, "gamma(-2.5) tight");
    ASSERT_NEAR(&ctx, ml_gamma_new(-3.5), 2.70088205852269114e-01,
                2.70088205852269114e-01 * 1.2e-15, "gamma(-3.5) tight");
    ASSERT_TRUE(&ctx, ml_isnan(ml_gamma_new(-1.0)), "gamma(-1) pole is NaN");
    ASSERT_TRUE(&ctx, ml_isnan(ml_gamma_new(-2.0)), "gamma(-2) pole is NaN");

    /* lgamma */
    ASSERT_TRUE(&ctx, ml_lgamma(1.0) == 0.0, "lgamma(1) == 0 exact");
    ASSERT_TRUE(&ctx, ml_lgamma(2.0) == 0.0, "lgamma(2) == 0 exact");
    ASSERT_NEAR(&ctx, ml_lgamma(5.0), ml_log(24.0), 1e-12, "lgamma(5) == log(24)");
    ASSERT_NEAR(&ctx, ml_lgamma(0.5), 5.72364942924700082e-01,
                5.72364942924700082e-01 * 1.2e-15, "lgamma(0.5) tight");
    ASSERT_NEAR(&ctx, ml_lgamma(1.5), -1.20782237635245218e-01,
                1.20782237635245218e-01 * 1.2e-15, "lgamma(1.5) tight");
    ASSERT_NEAR(&ctx, ml_lgamma(50.0), 1.44565743946344895e+02,
                1.44565743946344895e+02 * 1.2e-15, "lgamma(50) tight");
    ASSERT_NEAR(&ctx, ml_lgamma(100.0), 3.59134205369575398e+02,
                3.59134205369575398e+02 * 1.2e-15, "lgamma(100) tight");
    ASSERT_NEAR(&ctx, ml_lgamma(171.0), 7.06573062245787355e+02,
                7.06573062245787355e+02 * 1.2e-15, "lgamma(171) tight");
    ASSERT_NEAR(&ctx, ml_lgamma(-1.5), 8.60047015376480983e-01,
                8.60047015376480983e-01 * 1.2e-15, "lgamma(-1.5) tight");
    ASSERT_NEAR(&ctx, ml_lgamma(-2.5), -5.62437164976740539e-02,
                5.62437164976740539e-02 * 1.2e-15, "lgamma(-2.5) tight");
    ASSERT_TRUE(&ctx, ml_isinf(ml_lgamma(0.0)), "lgamma(0) is inf");
    ASSERT_TRUE(&ctx, ml_isnan(ml_lgamma(ml_make_nan())), "lgamma(NaN) is NaN");

    return ml_test_summary(&ctx);
}
"""


def fail(message: str) -> None:
    print("ERROR: " + message)
    sys.exit(1)


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


def patch_integral(v12: Path, force: bool) -> None:
    path = v12 / "src" / "integral.c"
    if not path.is_file():
        fail(f"Missing expected file: {path}")
    if MARKER in path.read_text(encoding="utf-8") and not force:
        print(f"  [skip] {path}: already at GAMMA_DD2")
        return
    write_text(path, NEW_INTEGRAL_C)


def patch_edge_test(v12: Path, force: bool) -> None:
    path = v12 / "tests" / "test_edge_integral.c"
    if not path.is_file():
        fail(f"Missing expected file: {path}")
    if TEST_MARKER in path.read_text(encoding="utf-8") and not force:
        print(f"  [skip] {path}: already at GAMMA_DD2_TEST")
        return
    write_text(path, NEW_EDGE_TEST_C)


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


def main() -> int:
    force = "--force" in sys.argv[1:]
    root, v12 = locate_v12a1()
    print("=========================================================")
    print("  MATHLIB v12A1: GAMMA DD FULL REWRITE (script 22)")
    print("=========================================================")
    print(f"  Root:  {root}")
    print(f"  v12A1: {v12}")
    print(f"  Force: {force}")
    print("---------------------------------------------------------")

    print("\n[1/2] integral.c — full rewrite (DD Lanczos)")
    patch_integral(v12, force)

    print("\n[2/2] test_edge_integral.c — oracle-aligned assertions")
    patch_edge_test(v12, force)

    archive_self(v12, force)

    print("\n---------------------------------------------------------")
    print("  Gamma DD rewrite applied.")
    print("")
    print("  Key differences from script 21:")
    print("    - DD low part reaches ml_exp: exp(L.hi) * (1 + L.lo)")
    print("    - true ~106-bit ml_log_dd (DD quotient + DD Horner)")
    print("    - DD Lanczos sum with DD term divisions")
    print("    - recurrence Gamma(x) = Gamma(x+1)/x for x < 0.5")
    print("    - whole-file rewrite: no dependence on 21's leftovers")
    print("")
    print("  Verify:")
    print("    cd v12A1")
    print("    cmake --build build")
    print("    ./build/test")
    print("    MATHLIB_EDGE_SANITIZERS=1 bash tests/run_edge_tests.sh")
    print("    ./build/oracle_check")
    print("")
    print("  Interpretation of the oracle afterwards:")
    print("    - integers 1..23 exact, large args single-digit ULP -> done")
    print("    - isolated 5-50 ULP residue at NON-integer points ->")
    print("      intrinsic error of the g=7,n=9 coefficient set;")
    print("      script 23 fits better coefficients with mpmath")
    print("      (or switches x >= 8 to DD Stirling).")
    print("=========================================================")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
