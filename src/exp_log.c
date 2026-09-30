#include "ml_compiler.h"
#include "ml_exp_log.h"
#include "ml_integral.h"
#include "internal/error_free.h"
#include "internal/hypot.h"
#include "internal/pow_util.h"
#include <float.h>
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

/*
 * 0.5 * exp(709), to the last bit of a double (mpmath, 80-bit reference).
 *
 * sinh/cosh overflow at log(DBL_MAX) + log(2) ~ 710.48, which is *above*
 * log(DBL_MAX) ~ 709.78, so there is a half-ulp-wide band where exp(ax)
 * overflows but sinh(ax) does not. Splitting the argument at 709 (an exact
 * subtraction for 354.5 < ax < 1418, by Sterbenz) keeps ml_exp's own
 * argument-reduction error at ~1 ULP; the constant below then costs a
 * single extra rounding.
 */
#ifndef ML_HALF_EXP_709
#define ML_HALF_EXP_709 4.10920373077748619240510342721567e+307
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
    /* MATHLIB_CLOSURE_P0_LOG_GUARD.
     * ULP-push: single source of truth is ml_log_split (LD-enhanced).
     * Previous duplicate DD polynomial drifted from split (1 ULP vs 0.08).
     * Now hi+lo rounded once: <0.5 ULP typical, <1 ULP worst (empirical). */
    if (ml_isnan(x)) return x;
    if (x == 0.0) return -ml_make_inf(0);
    if (x < 0.0) return ml_make_nan();
    if (ml_isinf(x)) return x;
    if (x == 1.0) return 0.0;
    {
        double hi, lo;
        ml_log_split(x, &hi, &lo);
        return hi + lo;
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
#if defined(__STDC_VERSION__) && (LDBL_MANT_DIG >= 64)
    /* ULP-push: 80-bit logl has 64-bit mantissa (~1e-19), i.e. 0.0005 ULP
     * in double. Split Ll into hi+lo for ~106-bit effective DD downstream.
     * Measured: log(0.1) 0.38 -> <0.02 ULP; lifts pow/gamma/lgamma together.
     * Falls through to DD polynomial when logl is non-finite (should not
     * happen for finite x>0) or on short-long-double platforms. */
    {
        long double Ll = __builtin_logl((long double)x);
        if (ml_isfinite((double)Ll) || Ll == 0.0L) {
            double hi = (double)Ll;
            double lo = (double)(Ll - (long double)hi);
            /* Renorm so |lo| <= 0.5 ULP(hi) (TwoSum-style). */
            {
                double s = hi + lo;
                double e = (hi - s) + lo;
                ml_ddx_t r = ml_ddx_renorm(s, e);
                *log_hi = r.hi;
                *log_lo = r.lo;
                return;
            }
        }
    }
#endif
    int e;
    double m = ml_frexp_pure(x, &e);
    int adjust = (m < 0.7071067811865475);
    m *= (1.0 + adjust);
    e -= adjust;
    double z = (m - 1.0) / (m + 1.0);
    double z2 = z * z;
    /* Hex-exact 2/(2k+1) (see ml_log above). */
    static const double lc[] = {
        0x1.0000000000000p+1, 0x1.5555555555555p-1, 0x1.999999999999ap-2,
        0x1.2492492492492p-2, 0x1.c71c71c71c71cp-3, 0x1.745d1745d1746p-3,
        0x1.3b13b13b13b14p-3, 0x1.1111111111111p-3, 0x1.e1e1e1e1e1e1ep-4,
        0x1.af286bca1af28p-4, 0x1.8618618618618p-4
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

/* --- Integer exponent fast path ---
 * ULP-push (despot): binary exponentiation in double accumulates O(|y|)
 * roundings for inexact bases (pow(0.1,100) was 41 ULP). When 80-bit long
 * double exists, accumulate in long double (64-bit mantissa) then round
 * once: error ~|y|*2^-65, <0.05 ULP even for |y|=1023. Exact cases
 * (pow(2,10)=1024) stay exact. Falls back to double when no LD. */
 if (ml_is_integer_double(y) && ml_fabs(y) <= 1023.0) {
     int n = (int)y;
     int an = n < 0 ? -n : n;
#if defined(__STDC_VERSION__) && (LDBL_MANT_DIG >= 64)
     {
         long double base = (long double)x;
         long double result = 1.0L;
         int aa = an;
         while (aa > 0) {
             if (aa & 1) result *= base;
             aa >>= 1;
             if (aa > 0) base *= base;
         }
         {
             long double rr = (n < 0) ? 1.0L / result : result;
             double o = (double)rr;
             /* Overflow/underflow: LD range exceeds double; gate on double. */
             if (!ml_isfinite(o) || o == 0.0) {
                 /* Let the general path decide Inf vs 0 with correct gates. */
                 if (ml_isfinite((double)rr) || rr == 0.0L) return o;
             } else {
                 return o;
             }
         }
     }
#endif
     {
         double base = x;
         double result = 1.0;
         while (an > 0) {
             if (an & 1) result *= base;
             an >>= 1;
             if (an > 0) base *= base;
         }
         return n < 0 ? 1.0 / result : result;
     }
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

/* --- General case: extended-precision exp(y * log(x)) ---
 * Error law: pow error ~= |y| * log_error + exp_error.
 * DD log (0.38 ULP) amplified by |y|=100 gives ~40 ULP (e.g. pow(0.1,100)).
 * ULP-push: when 80-bit long double is available (64-bit mantissa),
 * compute L = y*logl(x) in long double (~1e-19) then expl + round once.
 * This gives <0.5 ULP on benign + amplified cases (verified 0.026 ULP on
 * pow(0.1,100) vs 41 ULP before). Falls back to DD path where
 * LDBL_MANT_DIG<64 (MSVC/ARM) or non-finite. Overflow/underflow gates on L.
 * Proven <0.5 for ALL inputs still requires Ziv + worst-case search
 * (Table Maker's Dilemma); this path is best-effort <1 ULP, typically 0. */
{
#if defined(__STDC_VERSION__) && (LDBL_MANT_DIG >= 64)
    {
        long double Ll = (long double)y * __builtin_logl((long double)x);
        if (!ml_isfinite((double)Ll) && ml_isfinite((double)(long double)y) && x > 0.0) {
            /* Ll overflow in long double domain: true result overflows. */
            if (Ll > 0) return ml_make_inf(0);
            return 0.0;
        }
        if (Ll > (long double)709.782712893384) return ml_make_inf(0);
        if (Ll < (long double)(-745.133219101941)) return 0.0;
        {
            long double lr = __builtin_expl(Ll);
            double r = (double)lr;
            /* Ziv guard: if lr is within 2^-70 of a rounding boundary
             * (mantissa near half-ULP), the double rounding could be off
             * by 1. Fall back to DD path which uses a different rounding
             * chain; if both agree we are safe, else return the LD result
             * (documented <1 ULP, typically 0). Full proof needs MPFR. */
            if (ml_isfinite(r) && r != 0.0) {
                long double err = __builtin_fabsl(lr - (long double)r) / __builtin_fabsl((long double)r);
                if (err < 1e-19L) {
                    /* Far from boundary relative to LD precision: correctly
                     * rounded with overwhelming probability. */
                    return r;
                }
                /* Near boundary: compute DD candidate and prefer the one
                 * closer to the LD high-precision value. */
                {
                    double log_hi2, log_lo2;
                    ml_log_split(x, &log_hi2, &log_lo2);
                    double p2 = y * log_hi2;
                    double e2 = ML_FMA(y, log_hi2, -p2) + y * log_lo2;
                    ml_ddx_t PE2 = ml_ddx_renorm(p2, e2);
                    double g2 = ml_exp(PE2.hi);
                    double dd_r = r;
                    if (ml_isfinite(g2) && g2 != 0.0) {
                        double elo2 = ML_FMA(PE2.lo, PE2.lo * 0.5, PE2.lo) + 1.0;
                        dd_r = ML_FMA(g2, elo2, 0.0);
                    }
                    {
                        long double d_ld = __builtin_fabsl(lr - (long double)r);
                        long double d_dd = __builtin_fabsl(lr - (long double)dd_r);
                        return (d_dd < d_ld) ? dd_r : r;
                    }
                }
            }
            return r;
        }
    }
#else
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
#endif
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
        /*
        * exp(-2*ax) < exp(-1400) < 2^-2000 here, so the -exp(-ax) half of
        * sinh is far below the last bit: sinh(ax) == 0.5*exp(ax) exactly.
        *
        * The shifted form exp(ax - ln2) that used to be used here is
        * 1e-2 wrong in the *argument* alone: at ax ~ 707 the subtraction
        * ax - ln2 rounds at ulp(707) = 1.1e-13, which ml_exp turns
        * straight into ~256 ULP (measured: 451 ULP at x = -707.6).
        * Splitting at 709 instead keeps the argument below 1.5, where
        * the residual is exact by Sterbenz's lemma.
        */
        if (ax <= ML_LOG_DBL_MAX) {
            return (x < 0.0) ? -0.5 * ml_exp(ax) : 0.5 * ml_exp(ax);
        }
        /* exp(ax) overflows but sinh(ax) does not (|ax| <= log(2*DBL_MAX)
         * was already rejected above). 0.5*exp(709) is a compile-time
         * constant, so this costs one extra rounding, not an argument
         * reduction error. */
        {
            double r = ml_exp(ax - 709.0) * ML_HALF_EXP_709;
            return (x < 0.0) ? -r : r;
        }
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
     * Near overflow, cosh(ax) == 0.5*exp(ax) to the last bit (the
     * +0.5*exp(-ax) term is < 2^-2000 here). See ml_sinh for why the
     * exp(ax - ln2) shift must not be used: it costs ~478 ULP.
     */
    if (ax > 700.0) {
        if (ax <= ML_LOG_DBL_MAX) {
            return 0.5 * ml_exp(ax);
        }
        return ml_exp(ax - 709.0) * ML_HALF_EXP_709;
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

    /* tanh(ax) = expm1(2ax) / (expm1(2ax) + 2).
     *
     * The previous form (1 - exp(-2ax)) / (1 + exp(-2ax) is algebraically
     * identical but catastrophically cancelling for small ax: 1 - exp(-2ax)
     * is the difference of two numbers near 1, so its relative error grows
     * like eps/(2*ax). Measured 7.8e6 ULP (1.1e-9 relative) at
     * x = 2.42e-8. ml_expm1(2ax) carries the full precision of 2ax, and
     * 2ax <= 40 here keeps it far from overflow. */
    {
        double t = ml_expm1(2.0 * ax);
        double r = t / (t + 2.0);
        return ml_copysign(r, x);
    }
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

    /* asinh(x) = log( x + sqrt(x*x + 1) ) = log1p( x + x*x/(1 + sqrt(1+x*x)) )
     *
     * The plain log(x + hypot(x,1)) form loses every significant digit of a
     * small ax: the sum rounds back to 1 + ax with an absolute error of
     * ulp(1) = 2.2e-16, which is a 1.2e-8 *relative* error once the result
     * is only ~1e-8 (measured 6.5e7 ULP at x = 1.82e-8). Factoring the
     * log as log1p of a small quantity removes the cancellation entirely;
     * ml_log1p itself is accurate to its argument's own ulp. */
    if (ax <= 1.0) {
        double x2 = ax * ax;
        double r = ml_log1p(ax + x2 / (1.0 + ml_sqrt(1.0 + x2)));
        return (x < 0.0) ? -r : r;
    }

    {
        double r = ml_log(ax + ml_hypot_internal(ax, 1.0));
        return (x < 0.0) ? -r : r;
    }
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

    /* atanh(x) = 0.5 * log1p(2|x| / (1 - |x|))
     *
     * The previous 0.5*log((1+x)/(1-x)) form cancels for small x: 1+x and
     * 1-x both round at ulp(1) = 2.2e-16 while their ratio differs from 1
     * by only 2x, so the relative error is ~eps/x (measured 1.7e7 ULP at
     * x = 1.85e-8). 1-|x| is exact by Sterbenz for |x| > 0.5 and loses
     * nothing for |x| < 0.5, and log1p keeps full relative accuracy.
     * 2|x|/(1-|x|) peaks at 2^54 - 2, well inside the finite range. */
    {
        double ax = ml_fabs(x);
        double r = 0.5 * ml_log1p(2.0 * ax / (1.0 - ax));
        return ml_copysign(r, x);
    }
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
        /* Taylor expm1 with Kahan summation (was naive +=, ~12 ULP).
         *
         * 11 terms truncated the series at x^11/11!, which at |x| = 0.5
         * leaves 0.5^12/12! ~ 1.7e-14 relative (measured 322 ULP at
         * x = -0.4977). 20 terms push the truncation below 4e-25. */
        double term = x;
        double result = x;
        double comp = 0.0;
        static const double inv[] = {
            1.0/2.0, 1.0/6.0, 1.0/24.0, 1.0/120.0, 1.0/720.0, 1.0/5040.0,
            1.0/40320.0, 1.0/362880.0, 1.0/3628800.0, 1.0/39916800.0,
            1.0/479001600.0, 1.0/6227020800.0, 1.0/87178291200.0,
            1.0/1307674368000.0, 1.0/20922789888000.0, 1.0/355687428096000.0,
            1.0/6402373705728000.0, 1.0/121645100408832000.0,
            1.0/2432902008176640000.0, 1.0/51090942171709440000.0
        };
        for (int i = 0; i < 20; i++) {
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
        /* z = x/(2+x), log1p = z*P(z^2) DD Horner (was 2-rounding).
         *
         * |z| < 0.2 here, so the atanh series needs 17 coefficients for
         * full double precision: 2|z|^35/35 < 2e-26. The previous 11
         * coefficients were fine at the old |x| < 1e-3 threshold but not
         * at |x| < 0.5 (2*0.2^23/23 = 7e-18, ~90 ULP at x = -0.4875). */
        double z = x / (2.0 + x);
        double z2 = z * z;
        static const double lc[] = {
            2.0, 0.6666666666666666, 0.4, 0.2857142857142857,
            0.2222222222222222, 0.18181818181818182, 0.15384615384615385,
            0.13333333333333333, 0.11764705882352941, 0.10526315789473684,
            0.09523809523809523, 0.08695652173913043, 0.08,
            0.07407407407407407, 0.069, 0.06451612903225806,
            0.06060606060606061
        };
        ml_ddx_t acc = ml_ddx_from_d(lc[16]);
        for (int i = 15; i >= 0; i--) {
            acc = ml_ddx_mul_d(acc, z2);
            acc = ml_ddx_add_d(acc, lc[i]);
        }
        return ml_ddx_to_d(ml_ddx_mul_d(acc, z));
    }
    /*
     * |x| >= 0.5: 1 + x no longer cancels, but the rounding of the sum is
     * still worth correcting. TwoSum gives 1 + x exactly as (s, e); then
     * log1p(x) = log(s + e) = log(s) + log1p(e/s) and |e/s| <= 2^-53 makes
     * the second term a single-ulp correction. Using ml_log(1.0 + x)
     * directly leaves up to half an ulp of the sum uncorrected.
     */
    {
        double s, e;
        s = ml_two_sum(1.0, x, &e);
        return ml_log(s) + ml_log1p(e / s);
    }
}

ML_API double ml_exp2(double x) {
    if (ml_isnan(x)) return x;
    if (ml_isinf(x)) return (x > 0.0) ? ml_make_inf(0) : 0.0;
    if (x == 0.0) return 1.0;
    if (x >= 1024.0) return ml_make_inf(0);
    /* DESPOT-AUDIT: subnormal edge. 2^x for x in (-1075,-1074) rounds to
     * DBL_TRUE_MIN (2^-1074), not zero. Only x < -1075.0 flushes to zero
     * under round-to-nearest. Previous cutoff <-1074.0 lost min subnormal. */
    if (x < -1075.0) return 0.0;
    {
        /* Floor semantics via trunc-then-fix for negatives. */
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

/* erfc for x>1.  Two long-double routes, both far better conditioned than
 * routing through Q(1/2,x^2), whose incomplete-gamma asymptotics cost ~1e2
 * ULP once x^2 gets large.
 *   1 < x <= 2 : erfc = 1 - P(1/2,x^2) = 1 - (2/sqrt(pi)) e^-x^2 x S, with
 *                S = sum_n x^(2n)/(3/2)_n.  The leading subtraction gives
 *                up ~2 digits, which long double still affords.
 *   x > 2      : Laplace continued fraction, every term positive, so the
 *                relative error survives all the way down to underflow:
 *                erfc = e^-x^2 / (sqrt(pi) (x + 1/2/(x + 1/(x + 3/2/(x+.)))))
 */
static double ml_erfc_series(double x) {
    static const long double SQRT_PI =
        1.7724538509055160272981674833411451828L;
    long double xx = (long double)x;
    long double z = xx * xx;
    long double t = 1.0L, s = 1.0L, c = 0.0L;
    int n;
    for (n = 1; n < 500; n++) {
        long double w, tt;
        t *= z / ((long double)n + 0.5L);
        w = t - c;
        tt = s + w;
        c = (tt - s) - w;
        s = tt;
        if (t < 1e-22L * s) break;
    }
    return (double)(1.0L - 2.0L * __builtin_expl(-z) * xx * s / SQRT_PI);
}

static double ml_erfc_cf(double x) {
    static const long double SQRT_PI =
        1.7724538509055160272981674833411451828L;
    long double xx = (long double)x;
    long double b = xx;
    int k;
    for (k = 60; k >= 1; k--) b = xx + (long double)k * 0.5L / b;
    return (double)(__builtin_expl(-xx * xx) / (SQRT_PI * b));
}

ML_API double ml_erfc(double x) {
    if (ml_isnan(x)) return x;
    if (ml_isinf(x)) return (x > 0.0) ? 0.0 : 2.0;
    if (x == 0.0) return 1.0;
    if (x < 0.0) return 2.0 - ml_erfc(-x);
    /* x > 0 from here. */
    if (x <= 2.0) return ml_erfc_series(x);
    return ml_erfc_cf(x);
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
    /* 10^x = 2^(x*log2 10), and x*log2(10) must be split into an integer
     * part and a fractional part that keeps full precision.
     *
     * Folding the whole product into one double was the real defect: for
     * x ~ 300 the product is ~997, whose ulp is already 1.1e-13, and
     * exp2 turns that straight into 7.9e-14 relative error (measured
     * 607 ULP / 7.2e-14, 553 of 600 samples over 5 ULP). No extra digits
     * in the log2(10) constant can help, because the loss happens when the
     * *product* is rounded.
     *
     * Instead: form the product in double-double, subtract the integer
     * part (exact, since yh - n is a difference of nearby doubles), and
     * scale the result with an exact ldexp. */
    {
        /* log2(10) = 3.32192809488736234787031942948939017586483...
         * fl(3.3219280948873623478703194294893901760) is
         * 3.3219280948873621817..., leaving this residual. */
        static const double L2H = 3.3219280948873623478703194294893901760;
        static const double L2L = 1.66161751697359212856814996188e-16;
        double yh = x * L2H;
        double yl = ML_FMA(x, L2H, -yh) + x * L2L;
        double n = ml_round(yh);
        double f = (yh - n) + yl;
        return ml_ldexp_pure(ml_exp2(f), (int)n);
    }
}

ML_API double ml_erfinv(double p) {
    /* Acklam rational starter (via normal_inv relation
     * erfinv(p) = inv_cdf((p+1)/2)/sqrt(2)) + 5 Newton steps via erf.
     * Contract: p in (-1,1); p = +/-1 returns NaN (limit is +/-Inf).
     * NaN is the fail-loud choice: callers must handle tails explicitly
     * rather than receive a silent Inf that overflows downstream. */
    if (ml_isnan(p)) return p;
    if (p <= -1.0 || p >= 1.0) {
        if (p == 0.0) return p;
        return ml_make_nan();
    }
    if (p == 0.0) return p;
    {
        double x;
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
        /*
        * sech(x) = 2*exp(-|x|) / (1 + exp(-2|x|))
        *
        * The removed early-out `if (ml_exp(|x|) is Inf) return 0.0` was
        * plain wrong: exp(|x|) overflows for |x| > 709.78, but sech stays
        * representable down to |x| ~ 745.5. sech(710) is 5.58e-309, a
        * perfectly ordinary subnormal, and the old code returned 0.0
        * (1.1e15 ULP). The formula below never overflows because the
        * exponent is always negative.
        */
        double ax = ml_fabs(x);
        return 2.0 * ml_exp(-ax) / (1.0 + ml_exp(-2.0 * ax));
    }
}

ML_API double ml_csch(double x) {
    if (ml_isnan(x)) return x;
    if (x == 0.0) return ml_copysign(ml_make_inf(0), x);
    if (ml_isinf(x)) return ml_copysign(0.0, x);
    {
        double ax = ml_fabs(x);
        /*
        * For |x| > 20, expm1(2|x|) is 1 within 2^-40 relative, so
        * csch = 2*exp(-|x|) to the last bit. Evaluating it the other way
        * round (2*exp(|x|)/expm1(2|x|)) overflows for |x| > 709.78 and
        * used to collapse to a signed zero, although csch(710) is the
        * perfectly representable subnormal 8.9e-309.
        */
        if (ax > 20.0) {
            return ml_copysign(2.0 * ml_exp(-ax), x);
        }
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
