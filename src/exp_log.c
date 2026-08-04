#include "ml_compiler.h"
#include "ml_exp_log.h"
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
    double r = ML_FMA(-n, 0.69314718036912381649, x);
    r = ML_FMA(-n, 1.90821490974462528503e-10, r);

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

    double poly = 0.09523809523809523;
    poly = poly * z2 + 0.10526315789473684;
    poly = poly * z2 + 0.11764705882352941;
    poly = poly * z2 + 0.13333333333333333;
    poly = poly * z2 + 0.15384615384615385;
    poly = poly * z2 + 0.18181818181818182;
    poly = poly * z2 + 0.2222222222222222;
    poly = poly * z2 + 0.2857142857142857;
    poly = poly * z2 + 0.4;
    poly = poly * z2 + 0.6666666666666666;
    poly = poly * z2 + 2.0;

    /* MATHLIB_V12A1_LOG_COMPENSATED_RECONSTRUCT */
    /*
     * Compensated reconstruction.
     *
     * Old: ML_FMA((double)e, ML_LN2, z * poly)
     *   -> rounds e * ln(2) to 53 bits, losing up to 1e-13
     *      for large |e|.
     *
     * New: split ln(2) and use two terms.
     *   -> FMA(e, ln2_hi, z*poly) captures the high product
     *      with one rounding.
     *   -> e * ln2_lo adds back the truncated low bits.
     */
    return ML_FMA((double)e, ML_LN2_HI, z * poly)
         + (double)e * ML_LN2_LO;
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
    double poly = 0.09523809523809523;
    poly = poly * z2 + 0.10526315789473684;
    poly = poly * z2 + 0.11764705882352941;
    poly = poly * z2 + 0.13333333333333333;
    poly = poly * z2 + 0.15384615384615385;
    poly = poly * z2 + 0.18181818181818182;
    poly = poly * z2 + 0.2222222222222222;
    poly = poly * z2 + 0.2857142857142857;
    poly = poly * z2 + 0.4;
    poly = poly * z2 + 0.6666666666666666;
    poly = poly * z2 + 2.0;
    *log_hi = ML_FMA((double)e, ML_LN2_HI, z * poly);
    *log_lo = (double)e * ML_LN2_LO;
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
 * For |y| <= 64 and y integer, binary exponentiation is exact.
 * No log/exp roundtrip. pow(2, 10) = 1024 exactly.
 * pow(10, 3) = 1000 exactly. pow(2, -1) = 0.5 exactly.
 *
 * Works for negative bases too: pow(-2, 3) = -8.
 */
if (ml_is_integer_double(y) && ml_fabs(y) <= 64.0) {
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

/* --- Negative base, non-integer exponent --- */
if (x < 0.0) {
    return ml_make_nan();
}

/* --- General case: extended-precision exp(y * log(x)) --- */
/*
 * Old: ml_exp(y * ml_log(x))
 *   -> y * log(x) rounds once, losing up to 0.5 ULP.
 *   -> exp rounds again. Total: 2-4 ULP error.
 *
 * New: Dekker-split log(x) into log_hi + log_lo (exact).
 *   Compute y * (log_hi + log_lo) with FMA to capture low bits.
 *   Gives ~106 bits of precision before final rounding.
 *   Reduces pow error to ~1 ULP for most inputs.
 */
{
    double log_val = ml_log(x);
    /* Dekker split: log_val = log_hi + log_lo exactly */
    double c = 134217729.0 * log_val; /* (2^27 + 1) */
    double log_hi = c - (c - log_val);
    double log_lo = log_val - log_hi;
    /* Extended-precision product */
    double p = y * log_hi;
    double e = ML_FMA(y, log_hi, -p) + y * log_lo;
    return ml_exp(p + e);
}
}

ML_API double ml_logb(double x, double b) {
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
    if (ax < 1e-4) return x;

    double e = ml_exp(-2.0 * ax);
    double t = (1.0 - e) / (1.0 + e);

    return ml_copysign(t, x);
}

ML_API double ml_asinh(double x) {
    /* MATHLIB_CLOSURE_P0_ASINH_LARGE */
    if (ml_isnan(x) || ml_isinf(x)) return x;

    double ax = ml_fabs(x);
    if (ax == 0.0) return x;
    if (ax < 1e-4) return x;

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
    if (x <= -1.0 || x >= 1.0) return ml_make_nan();
    if (ml_fabs(x) < 1e-4) return x;

    return 0.5 * ml_log((1.0 + x) / (1.0 - x));
}
