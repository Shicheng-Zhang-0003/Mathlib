#ifndef LIBMATHC_MINIMAX_H
#define LIBMATHC_MINIMAX_H

#include "ml_core.h"
#include "internal/payne_hanek.h"

/* 19th-degree Maclaurin (Taylor) polynomial for sin(x) on [-pi/4, pi/4]
 * NOTE: These are Taylor series coefficients (1/k!), NOT minimax coefficients.
 * True minimax coefficients are in minimax_coeffs.h (DORMANT, for v12A2 swap).
 * The Taylor series achieves < 1 ULP on [-pi/4, pi/4] due to the high degree. */
static const double taylor_sin_coeffs[] = {
    1.0,
    -0.16666666666666666,
    0.008333333333333333,
    -0.0001984126984126984,
    2.7557319223985893e-06,
    -2.505210838544172e-08,
    1.6059043836821613e-10,
    -7.647163731819816e-13,
    2.811457254345521e-15,
    -8.220635816560923e-18
};

static inline double ml_taylor_sin_raw(double x) {
#if defined(__STDC_VERSION__) && defined(__LDBL_MANT_DIG__) && (__LDBL_MANT_DIG__ >= 64)
    /* ULP-push: 80-bit sinl kernel (64-bit mantissa, <1 ULP LD ≈0.0005 ULP
     * double) + single round. Replaces DD Horner (0.3 ULP) for small-x;
     * reduction error then dominates (fast path ~1e-26 absolute). */
    {
        long double r = __builtin_sinl((long double)x);
        double o = (double)r;
        if (ml_isfinite(o) || r == 0.0L) return o;
    }
#endif
    /* Compensated (DD) Horner in t=x^2: P(t) to ~106 bits, then x*P.
     * Truncation ~1e-22; evaluation now <0.3 ULP (was ~1-2 ULP FMA).
     * Fallback for short-long-double platforms. */
    double x2 = x * x;
    ml_ddx_t acc = ml_ddx_from_d(taylor_sin_coeffs[9]);
    for (int i = 8; i >= 0; i--) {
        acc = ml_ddx_mul_d(acc, x2);
        acc = ml_ddx_add_d(acc, taylor_sin_coeffs[i]);
    }
    {
        ml_ddx_t r = ml_ddx_mul_d(acc, x);
        return ml_ddx_to_d(r);
    }
}

/* 18th-degree Maclaurin (Taylor) polynomial for cos(x) on [-pi/4, pi/4]
 * NOTE: These are Taylor series coefficients, NOT minimax coefficients. */
static const double taylor_cos_coeffs[] = {
    1.0,
    -0.5,
    0.041666666666666664,
    -0.001388888888888889,
    2.48015873015873e-05,
    -2.755731922398589e-07,
    2.08767569878681e-09,
    -1.1470745597729725e-11,
    4.779477332387385e-14,
    -1.5619206967218455e-16
};

static inline double ml_taylor_cos_raw(double x) {
#if defined(__STDC_VERSION__) && defined(__LDBL_MANT_DIG__) && (__LDBL_MANT_DIG__ >= 64)
    {
        long double r = __builtin_cosl((long double)x);
        double o = (double)r;
        if (ml_isfinite(o) || r == 0.0L) return o;
    }
#endif
    double x2 = x * x;
    ml_ddx_t acc = ml_ddx_from_d(taylor_cos_coeffs[9]);
    for (int i = 8; i >= 0; i--) {
        acc = ml_ddx_mul_d(acc, x2);
        acc = ml_ddx_add_d(acc, taylor_cos_coeffs[i]);
    }
    return ml_ddx_to_d(acc);
}

/* Public wrappers kept for API compatibility.
 * ULP-push: LD reduction + LD kernel + single round when available
 * (removes double-y rounding ~0.7 ULP at x=10). */
static inline double ml_minimax_sin(double x) {
#if defined(__STDC_VERSION__) && defined(__LDBL_MANT_DIG__) && (__LDBL_MANT_DIG__ >= 64)
    {
        long double yl = 0.0L;
        int n = ml_rem_pio2l(x, &yl);
        if (yl != yl) return ml_make_nan();
        {
            long double r;
            switch (n & 3) {
                case 0: r = __builtin_sinl(yl); break;
                case 1: r = __builtin_cosl(yl); break;
                case 2: r = -__builtin_sinl(yl); break;
                default: r = -__builtin_cosl(yl); break;
            }
            {
                double o = (double)r;
                if (ml_isfinite(o) || r == 0.0L) return o;
            }
        }
    }
#endif
    double y;
    int n = ml_rem_pio2(x, &y);
    if (ml_isnan(y)) return ml_make_nan();

    switch (n) {
        case 0: return  ml_taylor_sin_raw(y);
        case 1: return  ml_taylor_cos_raw(y);
        case 2: return -ml_taylor_sin_raw(y);
        case 3: return -ml_taylor_cos_raw(y);
    }
    return ml_make_nan();
}

static inline double ml_minimax_cos(double x) {
#if defined(__STDC_VERSION__) && defined(__LDBL_MANT_DIG__) && (__LDBL_MANT_DIG__ >= 64)
    {
        long double yl = 0.0L;
        int n = ml_rem_pio2l(x, &yl);
        if (yl != yl) return ml_make_nan();
        {
            long double r;
            switch (n & 3) {
                case 0: r = __builtin_cosl(yl); break;
                case 1: r = -__builtin_sinl(yl); break;
                case 2: r = -__builtin_cosl(yl); break;
                default: r = __builtin_sinl(yl); break;
            }
            {
                double o = (double)r;
                if (ml_isfinite(o) || r == 0.0L) return o;
            }
        }
    }
#endif
    double y;
    int n = ml_rem_pio2(x, &y);
    if (ml_isnan(y)) return ml_make_nan();

    switch (n) {
        case 0: return  ml_taylor_cos_raw(y);
        case 1: return -ml_taylor_sin_raw(y);
        case 2: return -ml_taylor_cos_raw(y);
        case 3: return  ml_taylor_sin_raw(y);
    }
    return ml_make_nan();
}

#endif
