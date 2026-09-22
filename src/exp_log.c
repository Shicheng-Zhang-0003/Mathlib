#include "ml_compiler.h"
#include "ml_exp_log.h"
#include "ml_integral.h"
#include "internal/error_free.h"
#include "internal/hypot.h"
#include "internal/pow_util.h"
/* MATHLIB_CLOSURE_P2_P0_4_HYPERBOLIC_LIMITS */
#ifndef ML_LOG_DBL_MAX
#define ML_LOG_DBL_MAX 709.782712893384
#endif

#ifndef ML_LOG_HYP_OVERFLOW
#define ML_LOG_HYP_OVERFLOW (ML_LOG_DBL_MAX + ML_LN2)
#endif

/* MATHLIB_CLOSURE_P2_P0_3_EXP_LIMITS (ML_LOG_DBL_MAX defined above) */

#ifndef ML_LOG_UNDERFLOW
#define ML_LOG_UNDERFLOW (-745.133219101941)
#endif



/* MATHLIB_V12A1_LOG_COMPENSATED_RECONSTRUCT */
/*
 * Split ln(2) into high and low parts for compensated reconstruction.
 *
 * ML_LN2_HI has the low 26 bits of its significand zeroed.
 * ML_LN2_LO captures the remaining bits.
 * Together they represent ln(2) to ~106 bits of precision.
 *
 * These are the same values used by musl libc and match the
 * 2-term Cody-Waite split already used in ml_exp (script 03).
 */
#ifndef ML_LN2_HI
#define ML_LN2_HI 6.93147180369123816490e-01
#endif
#ifndef ML_LN2_LO
#define ML_LN2_LO 1.90821492927058500170e-10
#endif
/* v11S CLOSURE IP-4: overflow-safe hyperbolics */

ML_API double ml_exp(double x) {
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
    /* MATHLIB_V12A1_EXP_FMA_REDUCTION */
    /*
     * Error-free Cody-Waite reduction.
     *
     * The old code used two separate rounded subtractions:
     *   r = x - n * hi - n * lo
     *
     * Each multiply-subtract pair rounds twice, losing up to
     * 2 ULP of precision in the residual.
     *
     * FMA rounds once per step. The 2-term split of ln(2)
     * provides ~106 bits, which is sufficient for |n| <= 1024.
     */
    /* MATHLIB_V12A1_EXP_LN2_SPLIT_MACROS
     *
     * Use the canonical ln(2) split shared with ml_log().
     * The previous hardcoded low constant did not match ML_LN2_LO.
     */
    double r = ML_FMA(-n, ML_LN2_HI, x);
    r = ML_FMA(-n, ML_LN2_LO, r);

    static const double inv_fact[] = {
        1.0, 1.0, 0.5, 0.16666666666666666, 0.041666666666666664,
        0.008333333333333333, 0.001388888888888889, 0.0001984126984126984,
        2.48015873015873e-05, 2.7557319223985893e-06, 2.7557319223985888e-07,
        2.505210838544172e-08, 2.08767569878681e-09, 1.6059043836821613e-10,
        1.1470745597729725e-11, 7.647163731819816e-13, 4.779477332387385e-14,
        2.8114572543455206e-15, 1.5619206968586226e-16, 8.22063524662433e-18
    };

    /* DD Horner: P(r) to 106 bits, then exact ldexp. <0.5 ULP eval. */
    {
        ml_ddx_t acc = ml_ddx_from_d(inv_fact[19]);
        for (int i = 18; i >= 1; i--) {
            acc = ml_ddx_mul_d(acc, r);
            acc = ml_ddx_add_d(acc, inv_fact[i]);
        }
        acc = ml_ddx_mul_d(acc, r);
        acc = ml_ddx_add_d(acc, 1.0);
        {
            double p = ml_ddx_to_d(acc);
            int ni = (int)n;
            /* n==1024 at x~709.78: ldexp(p,1024)=Inf though true
             * exp(x)<=DBL_MAX. Split exponent to avoid spurious overflow. */
            if (ni > 1023) {
                return ml_ldexp_pure(p, 1023) * ml_ldexp_pure(1.0, ni - 1023);
            }
            return ml_ldexp_pure(p, ni);
        }
    }
}

ML_API double ml_log(double x) {
    /* MATHLIB_CLOSURE_P0_LOG_GUARD */
    if (ml_isnan(x)) return x;
    if (x == 0.0) return -ml_make_inf(0);
    if (x < 0.0) return ml_make_nan();
    if (ml_isinf(x)) return x;
    if (x == 1.0) return 0.0;

    int e;
    double m = ml_frexp_pure(x, &e);

    int adjust = (m < 0.7071067811865475);
    m *= (1.0 + adjust);
    e -= adjust;

    double z = (m - 1.0) / (m + 1.0);
    double z2 = z * z;

    /* DD Horner for atanh poly in t=z^2 (was 2-rounding mul+add). */
    static const double lc[] = {
        2.0, 0.6666666666666666, 0.4, 0.2857142857142857,
        0.2222222222222222, 0.18181818181818182, 0.15384615384615385,
        0.13333333333333333, 0.11764705882352941, 0.10526315789473684,
        0.09523809523809523
    };
    ml_ddx_t acc = ml_ddx_from_d(lc[10]);
    for (int i = 9; i >= 0; i--) {
        acc = ml_ddx_mul_d(acc, z2);
        acc = ml_ddx_add_d(acc, lc[i]);
    }
    {
        ml_ddx_t zp = ml_ddx_mul_d(acc, z);
        /* DD reconstruction: e*ln2 (DD) + zp. */
        double ed = (double)e;
        double ehi = ed * ML_LN2_HI;
        double elo = ML_FMA(ed, ML_LN2_HI, -ehi) + ed * ML_LN2_LO;
        ml_ddx_t eln2 = ml_ddx_renorm(ehi, elo);
        double s, e1;
        s = ml_two_sum(eln2.hi, zp.hi, &e1);
        e1 += eln2.lo + zp.lo;
        {
            ml_ddx_t r = ml_ddx_renorm(s, e1);
            return ml_ddx_to_d(r);
        }
    }
}
/* MATHLIB_V12A1_GAMMA_LOG_SPLIT */
/*
 * Double-double log: returns log(x) as log_hi + log_lo.
 *
 * log_hi is the main result (rounded to double).
 * log_lo captures the low bits of e*ln2.
 * Together they give ~106 bits of precision.
 *
 * This is used by ml_lgamma_positive and ml_gamma_new to avoid
 * the (z+0.5)*log(t) amplification error.
 */
ML_API void ml_log_split(double x, double *log_hi, double *log_lo) {
    if (ml_isnan(x) || x <= 0.0) {
        *log_hi = ml_make_nan();
        *log_lo = 0.0;
        return;
    }
    if (ml_isinf(x)) {
        *log_hi = x;
        *log_lo = 0.0;
        return;
    }
    if (x == 1.0) {
        *log_hi = 0.0;
        *log_lo = 0.0;
        return;
    }
    int e;
    double m = ml_frexp_pure(x, &e);
    int adjust = (m < 0.7071067811865475);
    m *= (1.0 + adjust);
    e -= adjust;
    double z = (m - 1.0) / (m + 1.0);
    double z2 = z * z;
    static const double lc[] = {
        2.0, 0.6666666666666666, 0.4, 0.2857142857142857,
        0.2222222222222222, 0.18181818181818182, 0.15384615384615385,
        0.13333333333333333, 0.11764705882352941, 0.10526315789473684,
        0.09523809523809523
    };
    ml_ddx_t acc = ml_ddx_from_d(lc[10]);
    for (int i = 9; i >= 0; i--) {
        acc = ml_ddx_mul_d(acc, z2);
        acc = ml_ddx_add_d(acc, lc[i]);
    }
    {
        ml_ddx_t zp = ml_ddx_mul_d(acc, z);
        double ed = (double)e;
        double ehi = ed * ML_LN2_HI;
        double elo = ML_FMA(ed, ML_LN2_HI, -ehi) + ed * ML_LN2_LO;
        ml_ddx_t eln2 = ml_ddx_renorm(ehi, elo);
        double s, e1;
        s = ml_two_sum(eln2.hi, zp.hi, &e1);
        e1 += eln2.lo + zp.lo;
        {
            ml_ddx_t r = ml_ddx_renorm(s, e1);
            *log_hi = r.hi;
            *log_lo = r.lo;
        }
    }
}


ML_API double ml_pow(double x, double y) {
/* MATHLIB_V12A1_POW_EXTENDED */

/* --- Special cases (unchanged from v11S) --- */
if (ml_isnan(y)) {
    if (x == 1.0) return 1.0;
    return ml_make_nan();
}
if (y == 0.0) return 1.0;
if (ml_isnan(x)) return ml_make_nan();
if (x == 1.0) return 1.0;

if (x == 0.0) {
    if (ml_isinf(y)) {
        return (y > 0.0) ? 0.0 : ml_make_inf(0);
    }
    if (y > 0.0) {
        if (ml_signbit(x) && ml_is_odd_integer_double(y)) {
            return ml_copysign(0.0, -1.0);
        }
        return 0.0;
    }
    if (ml_signbit(x) && ml_is_odd_integer_double(y)) {
        return -ml_make_inf(0);
    }
    return ml_make_inf(0);
}

if (ml_isinf(y)) {
    double ax = ml_fabs(x);
    if (ax == 1.0) return 1.0;
    if (y > 0.0) return (ax > 1.0) ? ml_make_inf(0) : 0.0;
    return (ax > 1.0) ? 0.0 : ml_make_inf(0);
}

if (ml_isinf(x)) {
    if (x > 0.0) {
        return (y > 0.0) ? ml_make_inf(0) : 0.0;
    }
    if (!ml_is_integer_double(y)) return ml_make_nan();
    if (y > 0.0) {
        return ml_is_odd_integer_double(y)
             ? -ml_make_inf(0) : ml_make_inf(0);
    }
    return ml_is_odd_integer_double(y)
         ? ml_copysign(0.0, -1.0) : 0.0;
}

/* --- Integer exponent fast path --- */
 /*
  * For |y| <= 1023 and y integer, binary exponentiation is exact.
  * No log/exp roundtrip. pow(2, 10) = 1024 exactly.
  * pow(10, 3) = 1000 exactly. pow(2, -1) = 0.5 exactly.
  *
  * Works for negative bases too: pow(-2, 3) = -8.
  */
 if (ml_is_integer_double(y) && ml_fabs(y) <= 1023.0) {
     int n = (int)y;
     int an = n < 0 ? -n : n;
     double base = x;
     double result = 1.0;
     while (an > 0) {
         if (an & 1) result *= base;
         an >>= 1;
         if (an > 0) base *= base;
     }
     return n < 0 ? 1.0 / result : result;
 }

/* --- Negative base, integer exponent (any magnitude) --- */
/* Integer-valued y with negative x is defined: sign = (-1)^y.
 * The |y|<=1023 fast path above is exact; for larger |y| use
 * exp(y*log|x|) for magnitude then apply odd/even sign.
 * This fixes pow(-2,2000)=+Inf and pow(-1,1e308)=1, which
 * previously fell through to NaN. */
if (x < 0.0) {
    if (ml_is_integer_double(y)) {
        double mag = ml_pow(-x, y);
        if (ml_isnan(mag)) return ml_make_nan();
        if (ml_is_odd_integer_double(y)) mag = -mag;
        /* Preserve signed-zero semantics: (-0)^odd etc. handled
         * by the x==0 path above; here x<0 and finite nonzero. */
        return mag;
    }
    return ml_make_nan();
}

/* --- General case: DD exp(y * log(x)) for <1 ULP ---
 * log via DD split (106 bits), product y*log as DD (p,e),
 * then exp_dd (exp(hi)*exp(lo) 2nd-order) instead of rounding
 * p+e to double first (which cost 0.5 ULP). */
{
    double log_hi, log_lo;
    ml_log_split(x, &log_hi, &log_lo);
    double p = y * log_hi;
    double e = ML_FMA(y, log_hi, -p) + y * log_lo;
    {
        ml_ddx_t PE = ml_ddx_renorm(p, e);
        if (PE.hi > 709.782712893384) return ml_make_inf(0);
        if (PE.hi < -745.133219101941) return 0.0;
        {
            double g = ml_exp(PE.hi);
            if (!ml_isfinite(g) || g == 0.0) return g;
            {
                double elo = ML_FMA(PE.lo, PE.lo * 0.5, PE.lo) + 1.0;
                return ML_FMA(g, elo, 0.0);
            }
        }
    }
}
}

ML_API double ml_logb(double x, double b) {
    /* log_b(x) undefined for b==1 (division by log(1)=0), b<=0,
     * x<=0 (except IEEE edge handling already in ml_log). */
    if (ml_isnan(x) || ml_isnan(b)) return ml_make_nan();
    if (b == 1.0) return ml_make_nan();
    if (b <= 0.0) return ml_make_nan();
    return ml_log(x) / ml_log(b);
}

ML_API double ml_sinh(double x) {
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
}
ML_API double ml_cosh(double x) {
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

ML_API double ml_tanh(double x) {
    if (ml_isnan(x)) return x;
    if (ml_isinf(x)) return ml_copysign(1.0, x);

    double ax = ml_fabs(x);

    if (ax == 0.0) return x;
    if (ax > 20.0) return ml_copysign(1.0, x);
    if (ax < 1.5e-8) return x;

    double e = ml_exp(-2.0 * ax);
    double t = (1.0 - e) / (1.0 + e);

    return ml_copysign(t, x);
}

ML_API double ml_asinh(double x) {
    /* MATHLIB_CLOSURE_P0_ASINH_LARGE */
    if (ml_isnan(x) || ml_isinf(x)) return x;

    double ax = ml_fabs(x);
    if (ax == 0.0) return x;
    if (ax < 1.5e-8) return x;

    if (ax > 1e150) {
        double r = ml_log(2.0) + ml_log(ax);
        return (x < 0.0) ? -r : r;
    }

    double r = ml_log(ax + ml_hypot_internal(ax, 1.0));
    return (x < 0.0) ? -r : r;
}

ML_API double ml_acosh(double x) {
    if (ml_isnan(x)) return x;
    if (x < 1.0) return ml_make_nan();
    if (x == 1.0) return 0.0;
    if (ml_isinf(x)) return x;

    if (x > 1e150) {
        return ml_log(2.0) + ml_log(x);
    }

    return ml_log(x + ml_sqrt((x - 1.0) * (x + 1.0)));
}

ML_API double ml_atanh(double x) {
    if (ml_isnan(x)) return x;
    if (x == 1.0) return ml_make_inf(0);
    if (x == -1.0) return ml_make_inf(1);
    if (x < -1.0 || x > 1.0) return ml_make_nan();
    if (ml_fabs(x) < 1.5e-8) return x;

    return 0.5 * ml_log((1.0 + x) / (1.0 - x));
}

/* ---- libm completion: expm1/log1p/exp2/log2/log10/cbrt/erf ---- */

ML_API double ml_expm1(double x) {
    if (ml_isnan(x)) return x;
    if (ml_isinf(x)) return (x > 0.0) ? ml_make_inf(0) : -1.0;
    if (x == 0.0) return x;
    double ax = ml_fabs(x);
    /* Return-x threshold is QUADRATIC (x^2/2): |x|<eps/2 for <0.5 ULP.
     * (1.5e-8 is the CUBIC threshold for sinh/tanh, NOT valid here —
     *  it cost 2.4M ULP at x=1e-9.) */
    if (ax < 1.11e-16) return x;
    if (ax < 0.5) {
        /* Taylor expm1 with Kahan summation (was naive +=, ~12 ULP). */
        double x2 = x * x;
        (void)x2;
        double term = x;
        double result = x;
        double comp = 0.0;
        static const double inv[] = {
            0.5, 1.0/6.0, 1.0/24.0, 1.0/120.0, 1.0/720.0, 1.0/5040.0,
            1.0/40320.0, 1.0/362880.0, 1.0/3628800.0, 1.0/39916800.0,
            1.0/479001600.0
        };
        for (int i = 0; i < 11; i++) {
            term *= x;
            {
                double w = term * inv[i] - comp;
                double t = result + w;
                comp = (t - result) - w;
                result = t;
            }
        }
        return result;
    }
    {
        double e = ml_exp(x);
        if (ml_isinf(e)) return e;
        return e - 1.0;
    }
}

ML_API double ml_log1p(double x) {
    if (ml_isnan(x)) return x;
    if (x == 0.0) return x;
    if (x < -1.0) return ml_make_nan();
    if (x == -1.0) return -ml_make_inf(0);
    if (ml_isinf(x)) return x;
    double ax = ml_fabs(x);
    /* Quadratic error x^2/2: same eps/2 threshold as expm1. */
    if (ax < 1.11e-16) return x;
    if (ax < 0.5) {
        /* z = x/(2+x), log1p = z*P(z^2) DD Horner (was 2-rounding). */
        double z = x / (2.0 + x);
        double z2 = z * z;
        static const double lc[] = {
            2.0, 0.6666666666666666, 0.4, 0.2857142857142857,
            0.2222222222222222, 0.18181818181818182, 0.15384615384615385,
            0.13333333333333333, 0.11764705882352941, 0.10526315789473684,
            0.09523809523809523
        };
        ml_ddx_t acc = ml_ddx_from_d(lc[10]);
        for (int i = 9; i >= 0; i--) {
            acc = ml_ddx_mul_d(acc, z2);
            acc = ml_ddx_add_d(acc, lc[i]);
        }
        return ml_ddx_to_d(ml_ddx_mul_d(acc, z));
    }
    return ml_log(1.0 + x);
}

ML_API double ml_exp2(double x) {
    if (ml_isnan(x)) return x;
    if (ml_isinf(x)) return (x > 0.0) ? ml_make_inf(0) : 0.0;
    if (x == 0.0) return 1.0;
    if (x >= 1024.0) return ml_make_inf(0);
    if (x < -1074.0) return 0.0;
    {
        double n = ml_round(x);
        /* Use floor semantics via round-then-fix for negatives. */
        double fl = (x >= 0.0) ? ml_trunc(x) : (ml_trunc(x) - ((x == ml_trunc(x)) ? 0.0 : 1.0));
        double f = x - fl;
        double p;
        if (f == 0.0) {
            p = 1.0;
        } else {
            /* DD product f*ln2 (was 1 rounding), then DD exp. */
            double ph = f * ML_LN2_HI;
            double pl = ML_FMA(f, ML_LN2_HI, -ph) + f * ML_LN2_LO;
            ml_ddx_t PE = ml_ddx_renorm(ph, pl);
            double g = ml_exp(PE.hi);
            if (!ml_isfinite(g)) {
                p = ml_exp(f * ML_LN2);
            } else {
                double elo = ML_FMA(PE.lo, PE.lo * 0.5, PE.lo) + 1.0;
                p = ML_FMA(g, elo, 0.0);
            }
        }
        (void)n;
        return ml_ldexp_pure(p, (int)fl);
    }
}

ML_API double ml_log2(double x) {
    if (ml_isnan(x)) return x;
    if (x <= 0.0) {
        if (x == 0.0) return -ml_make_inf(0);
        return ml_make_nan();
    }
    if (ml_isinf(x)) return x;
    if (x == 1.0) return 0.0;
    {
        double hi, lo;
        ml_log_split(x, &hi, &lo);
        /* DD division by ln2: q=hi/L, q+=(FMA(-q,L,hi)+lo)/L. */
        double q = hi / ML_LN2;
        double r = ML_FMA(-q, ML_LN2, hi) + lo;
        q += r / ML_LN2;
        return q;
    }
}

ML_API double ml_log10(double x) {
    static const double LN10_HI = 2.3025850929940459;
    static const double LN10_LO = 2.1707562233822494e-16;
    if (ml_isnan(x)) return x;
    if (x <= 0.0) {
        if (x == 0.0) return -ml_make_inf(0);
        return ml_make_nan();
    }
    if (ml_isinf(x)) return x;
    if (x == 1.0) return 0.0;
    {
        double hi, lo;
        ml_log_split(x, &hi, &lo);
        double q = hi / LN10_HI;
        double r = ML_FMA(-q, LN10_HI, hi) + lo - q * LN10_LO;
        q += r / LN10_HI;
        return q;
    }
}

ML_API double ml_cbrt(double x) {
    /* Builtin + one Newton refinement for perfect cubes (27->3 exact).
     * Builtin alone is 1 ULP high on exact cubes due to runtime rounding. */
    if (ml_isnan(x) || x == 0.0 || ml_isinf(x)) return x;
    {
        double r = __builtin_cbrt(x);
        if (ml_isfinite(r) && r != 0.0) {
            double r2 = r * r;
            double q = x / r2;
            double nr = (2.0 * r + q) / 3.0;
            if (ml_isfinite(nr)) {
                /* Accept refinement only if it reduces |r^3-x|. */
                double e0 = ml_fabs(ML_FMA(r2, r, -x));
                double n2 = nr * nr;
                double e1 = ml_fabs(ML_FMA(n2, nr, -x));
                if (!(e1 > e0)) r = nr;
            }
        }
        return r;
    }
}

/* erf/erfc: Taylor(|x|<=1, Kahan, ~1 ULP) + erf(1)+quad integral(1<|x|<=4)
 * + asymptotic(|x|>4). No external tables; no change to existing kernels. */
static double ml_erf_taylor(double x) {
    /* erf = 2/sqrt(pi) * sum (-1)^n x^{2n+1}/(n!(2n+1)), Kahan. */
    static const double C = 1.12837916709551257390; /* 2/sqrt(pi) */
    double x2 = x * x;
    double term = x;
    double sum = x;
    double comp = 0.0;
    for (int n = 1; n < 60; n++) {
        term *= -x2 / (double)n;
        {
            double w = term / (double)(2 * n + 1) - comp;
            double t = sum + w;
            comp = (t - sum) - w;
            sum = t;
        }
        if (ml_fabs(term / (double)(2 * n + 1)) < 1e-19 * ml_fabs(sum) + 1e-300) break;
    }
    return C * sum;
}

ML_API double ml_erfc(double x) {
    if (ml_isnan(x)) return x;
    if (ml_isinf(x)) return (x > 0.0) ? 0.0 : 2.0;
    if (x == 0.0) return 1.0;
    if (x < 0.0) return 2.0 - ml_erfc(-x);
    /* x > 0 from here. */
    if (x <= 1.0) {
        return 1.0 - ml_erf_taylor(x);
    }
    /* x > 1: erfc(x) = Q(1/2,x^2) via incomplete gamma (few ULP).
     * Replaces 512-pt Simpson (~500 ULP from 512 exp roundings). */
    return ml_gamma_q(0.5, x * x);
}

ML_API double ml_erf(double x) {
    if (ml_isnan(x) || x == 0.0 || ml_isinf(x)) {
        if (ml_isnan(x)) return x;
        if (x == 0.0) return x;
        return (x > 0.0) ? 1.0 : -1.0;
    }
    {
        double ax = ml_fabs(x);
        if (ax <= 1.0) return ml_erf_taylor(x);
        /* |x|>1: erf = sign*P(1/2,x^2). Few ULP via gamma CF/series. */
        {
            double p = ml_gamma_p(0.5, ax * ax);
            return (x > 0.0) ? p : -p;
        }
    }
}

static double ml_lambert_iter(double x, double w) {
    /* One Fritsch/Halley step: w -= 2(we^w-x)/(...). Caller loops. */
    double ew = ml_exp(w);
    double f = w * ew - x;
    double denom = ew * (w + 1.0) - (w + 2.0) * f / (2.0 * w + 2.0);
    if (denom == 0.0 || !ml_isfinite(denom)) return ml_make_nan();
    return w - f / denom;
}

ML_API double ml_lambert_w0(double x) {
    static const double INV_E = 0.36787944117144232160;
    if (ml_isnan(x)) return x;
    if (x == 0.0) return x;
    if (ml_isinf(x)) return (x > 0.0) ? x : ml_make_nan();
    if (x < -INV_E) return ml_make_nan();
    if (x == -INV_E) return -1.0;
    {
        double w;
        if (x < 3.0) {
            /* log1p-based start: w≈ln(1+x) within 0.2 of root for x<3,
             * avoids L1=0 singularity at x=1 in L1-L2 form. */
            w = ml_log1p(x);
            if (!ml_isfinite(w)) w = 0.5;
        } else {
            double L1 = ml_log(x), L2 = ml_log(L1);
            w = L1 - L2 + L2 / L1;
        }
        for (int i = 0; i < 20; i++) {
            double nw = ml_lambert_iter(x, w);
            if (!ml_isfinite(nw)) return ml_make_nan();
            if (ml_fabs(nw - w) < 1e-16 * (1.0 + ml_fabs(nw))) return nw;
            w = nw;
        }
        return w;
    }
}

ML_API double ml_lambert_wm1(double x) {
    static const double INV_E = 0.36787944117144232160;
    if (ml_isnan(x) || x >= 0.0) {
        if (ml_isnan(x)) return x;
        return ml_make_nan();
    }
    if (x < -INV_E) return ml_make_nan();
    if (x == -INV_E) return -1.0;
    {
        double L1 = ml_log(-x), L2 = ml_log(-L1);
        double w = L1 - L2 + L2 / L1;
        for (int i = 0; i < 30; i++) {
            double nw = ml_lambert_iter(x, w);
            if (!ml_isfinite(nw)) return ml_make_nan();
            if (ml_fabs(nw - w) < 1e-16 * (1.0 + ml_fabs(nw))) return nw;
            w = nw;
        }
        return w;
    }
}

ML_API double ml_exp10(double x) {
    if (ml_isnan(x)) return x;
    if (ml_isinf(x)) return (x > 0.0) ? ml_make_inf(0) : 0.0;
    if (x == 0.0) return 1.0;
    if (x >= 309.0) return ml_make_inf(0);
    /* 10^-323 ~ 9.9e-324 is a finite subnormal (min ~4.9e-324 = 10^-323.3);
     * old -323.0 cutoff flushed it to 0. Route via exp2/ldexp instead. */
    if (x < -324.0) return 0.0;
    /* 10^x = 2^{x log2 10}; exp2 does exact integer ldexp. */
    {
        static const double L2 = 3.32192809488736234787;
        double y = x * L2;
        /* One FMA refinement of y = x*L2 with LO of log2(10):
         * LO = 3.4e-17 (exact hex split). */
        y += x * 3.41765073873262459376e-17;
        return ml_exp2(y);
    }
}

ML_API double ml_erfinv(double p) {
    /* Winitzki/Strecok initial + 3 Halley steps via erf. */
    if (ml_isnan(p)) return p;
    if (p <= -1.0 || p >= 1.0) {
        if (p == 0.0) return p;
        return ml_make_nan();
    }
    if (p == 0.0) return p;
    {
        double w = -ml_log((1.0 - p) * (1.0 + p));
        double x;
        (void)w;
        /* Acklam-style central starter via normal_inv relation:
         * erfinv(p) = inv_cdf((p+1)/2)/sqrt(2). Reuse rational below. */
        {
            double q = (p + 1.0) * 0.5;
            double s = ml_sqrt(2.0);
            /* Use ml_normal_inv equivalent without mu/sigma dependency:
             * inline Acklam to avoid header cycle (statistics includes us). */
            double z;
            if (q < 0.02425) {
                double qq = ml_sqrt(-2.0 * ml_log(q));
                z = (((((-7.784894002430293e-03 * qq - 0.3223964580411365) * qq - 2.400758277161838) * qq - 2.549732539343734) * qq + 4.374664141464968) * qq + 2.938163982698783) / ((((7.784695709041462e-03 * qq + 0.3224671290700398) * qq + 2.445134137142996) * qq + 3.754408661907416) * qq + 1.0);
            } else if (q > 0.97575) {
                double qq = ml_sqrt(-2.0 * ml_log(1.0 - q));
                z = -(((((-7.784894002430293e-03 * qq - 0.3223964580411365) * qq - 2.400758277161838) * qq - 2.549732539343734) * qq + 4.374664141464968) * qq + 2.938163982698783) / ((((7.784695709041462e-03 * qq + 0.3224671290700398) * qq + 2.445134137142996) * qq + 3.754408661907416) * qq + 1.0);
            } else {
                double qq = q - 0.5, rr = qq * qq;
                z = (((((-39.69683028665376 * rr + 220.9460984245205) * rr - 275.9285104469687) * rr + 138.3577518672690) * rr - 30.66479806614716) * rr + 2.506628277459239) * qq / (((((-54.47609879822406 * rr + 161.5858368580409) * rr - 155.6989798598866) * rr + 66.80131188771972) * rr - 13.28068155288572) * rr + 1.0);
            }
            x = z / s;
        }
        for (int i = 0; i < 5; i++) {
            double e = ml_erf(x) - p;
            double d = ml_exp(-x * x) * 1.12837916709551257390;
            if (d == 0.0 || !ml_isfinite(d)) break;
            {
                double dx = e / d;
                x -= dx;
                if (ml_fabs(dx) < 1e-16 * (1.0 + ml_fabs(x))) break;
            }
        }
        return x;
    }
}

ML_API double ml_sech(double x) {
    if (ml_isnan(x)) return x;
    if (ml_isinf(x)) return 0.0;
    {
        double e = ml_exp(ml_fabs(x));
        if (ml_isinf(e)) return 0.0;
        return 2.0 * ml_exp(-ml_fabs(x)) / (1.0 + ml_exp(-2.0 * ml_fabs(x)));
    }
}

ML_API double ml_csch(double x) {
    if (ml_isnan(x)) return x;
    if (x == 0.0) return ml_copysign(ml_make_inf(0), x);
    if (ml_isinf(x)) return ml_copysign(0.0, x);
    {
        double ax = ml_fabs(x);
        double e = ml_expm1(2.0 * ax);
        if (!ml_isfinite(e)) return ml_copysign(0.0, x);
        if (e == 0.0) {
            /* ax subnormal-small: 1-exp(-2ax) would round to 0.
             * csch(ax)~1/ax via expm1 path. */
            return ml_copysign(1.0 / ax - ax / 6.0, x);
        }
        /* 2*exp(-ax)/(1-exp(-2ax)) = 2*exp(ax)/expm1(2ax): uses the
         * accurate expm1 instead of the rounded 1-exp(-2ax). */
        double r = 2.0 * ml_exp(ax) / e;
        return ml_copysign(r, x);
    }
}

ML_API double ml_coth(double x) {
    if (ml_isnan(x)) return x;
    if (x == 0.0) return ml_copysign(ml_make_inf(0), x);
    if (ml_isinf(x)) return ml_copysign(1.0, x);
    {
        double ax = ml_fabs(x);
        if (ax < 1e-16) {
            /* (1+e)/(1-e) would divide by rounded 0; coth~1/ax. */
            return ml_copysign(1.0 / ax, x);
        }
        double e = ml_exp(-2.0 * ax);
        /* Denominator via -expm1(-2ax) so tiny ax stays nonzero. */
        double denom = -ml_expm1(-2.0 * ax);
        if (denom == 0.0) return ml_copysign(ml_make_inf(0), x);
        double r = (1.0 + e) / denom;
        /* coth = (1+e)/(1-e); 1+e computed directly is exact here
         * since e in (0,1); denom above is the accurate part. */
        (void)e;
        return ml_copysign(r, x);
    }
}

ML_API double ml_asech(double x) {
    if (ml_isnan(x)) return x;
    if (x <= 0.0 || x > 1.0) return ml_make_nan();
    if (x == 1.0) return 0.0;
    {
        /* log((1+s)/x) = log(1+s)-log(x): avoids 2/x overflow for
         * x<1e-308 (would give Inf instead of ~744). */
        double s = ml_sqrt(1.0 - x * x);
        return ml_log(1.0 + s) - ml_log(x);
    }
}

ML_API double ml_acsch(double x) {
    if (ml_isnan(x)) return x;
    if (x == 0.0) return ml_copysign(ml_make_inf(0), x);
    if (ml_isinf(x)) return ml_copysign(0.0, x);
    {
        double ax = ml_fabs(x);
        if (ax < 1e-150) {
            /* 1/(x*x) overflows below ~1e-154; use log(2/|x|) tail:
             * acsch(x)~sign*log(2/|x|) to <1 ULP for tiny x. */
            return ml_copysign(ml_log(2.0) - ml_log(ax), x);
        }
        return ml_log(1.0 / x + ml_sqrt(1.0 / (x * x) + 1.0));
    }
}

ML_API double ml_acoth(double x) {
    if (ml_isnan(x)) return x;
    if (x == 1.0) return ml_make_inf(0);
    if (x == -1.0) return ml_make_inf(1);
    if (ml_fabs(x) <= 1.0) return ml_make_nan();
    if (ml_isinf(x)) return ml_copysign(0.0, x);
    if (ml_fabs(x) > 2.0) {
        /* (x+1)/(x-1) rounds to 1 for |x|>>1 -> log=0 instead of ~1/x.
         * log1p(2/(x-1))/2 stays accurate. */
        return 0.5 * ml_log1p(2.0 / (x - 1.0));
    }
    return 0.5 * ml_log((x + 1.0) / (x - 1.0));
}
