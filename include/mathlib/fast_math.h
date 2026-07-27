#ifndef LIBMATHC_FAST_MATH_H
#define LIBMATHC_FAST_MATH_H

#include <stdint.h>
#include <string.h>
#include "ml_core.h"
#include "bitwise_fp.h"

/* ============================================================================
 * FAST MATH APPROXIMATIONS (v11S CLOSURE HARDENING)
 *
 * These functions remain approximate, but they are no longer allowed to:
 *   - invoke undefined behavior
 *   - silently mishandle subnormals
 *   - produce garbage on NaN / Inf / zero inputs
 *
 * Note:
 * Uses memcpy instead of union type-punning to guarantee ISO C99 compliance
 * and prevent strict-aliasing undefined behavior.
 * ========================================================================== */

/*
 * Fast inverse square root for 64-bit doubles.
 *
 * Contract:
 *   - NaN input -> NaN
 *   - negative finite input -> NaN
 *   - +Inf -> 0
 *   - zero -> signed Inf
 *   - positive finite -> approximate rsqrt
 */
ML_INLINE double ml_fast_rsqrt(double number) {
    /* v11S AUDIT IP-1: positive subnormal fallback */
    if (ml_isnan(number)) return number;
    if (number < 0.0) return ml_make_nan();
    if (ml_isinf(number)) return 0.0;
    if (number == 0.0) return ml_copysign(ml_make_inf(0), number);

    if (ml_is_subnormal(number)) {
        return 1.0 / ml_sqrt(number);
    }

    uint64_t bits;
    memcpy(&bits, &number, sizeof(uint64_t));

    bits = 0x5fe6ec85e7de30daULL - (bits >> 1);

    double y;
    memcpy(&y, &bits, sizeof(double));

    y = y * (1.5 - (number * 0.5 * y * y));
    return y * (1.5 - (number * 0.5 * y * y));
}

/*
 * Fast log2 using frexp-based normalization.
 *
 * This is now subnormal-aware because it routes through ml_frexp_pure()
 * instead of directly reading the raw exponent field.
 *
 * Contract:
 *   - NaN -> NaN
 *   - negative -> NaN
 *   - zero -> -Inf
 *   - +Inf -> +Inf
 *   - positive finite -> approximate log2(x)
 */
ML_INLINE double ml_fast_log2(double x) {
    if (ml_isnan(x)) return x;
    if (x < 0.0) return ml_make_nan();
    if (x == 0.0) return -ml_make_inf(0);
    if (ml_isinf(x)) return x;

    int e;
    double m = ml_frexp_pure(x, &e);

    /*
     * m is in [0.5, 1.0)
     * log2(x) = e + log2(m)
     *
     * Let f = 2m - 1, so f is in [0, 1)
     * log2(m) = -1 + log2(1+f)
     */
    double f = 2.0 * m - 1.0;

    /* 3rd-degree polynomial for log2(1+f) on [0, 1] */
    double p = f * (1.442695 + f * (-0.721347 + f * 0.278652));

    return (double)(e - 1) + p;
}

/*
 * Fast exp2 using reverse isomorphism + fractional polynomial.
 *
 * v11S CLOSURE FIX:
 * The old implementation cast x directly to long long without bounds,
 * which is undefined behavior for huge inputs.
 *
 * Contract:
 *   - NaN -> NaN
 *   - x >= 1024 -> +Inf
 *   - x < -1074 -> 0
 *   - otherwise -> approximate 2^x
 */
ML_INLINE double ml_fast_exp2(double x) {
    if (ml_isnan(x)) return x;
    if (x >= 1024.0) return ml_make_inf(0);
    if (x < -1074.0) return 0.0;

    /* Safe: x is now bounded. */
    long long xi = (long long)x;

    /* Floor behavior for negative non-integers */
    if (x < 0.0 && (double)xi != x) xi -= 1;

    double exp_int = (double)xi;
    double frac = x - exp_int;

    /* 5th-degree polynomial for 2^frac on [0, 1] */
    double p = frac * (
        ML_LN2 + frac * (
            0.2402265069591007 + frac * (
                0.05550410866482158 + frac * (
                    0.009618129107628477 + frac * 0.00133335581464249
                )
            )
        )
    );

    double mant_approx = 1.0 + p;

    long long biased = xi + 1023;

    /*
     * If biased exponent is non-positive, we are in subnormal / underflow
     * territory. Route through ml_ldexp_pure() so the core IEEE layer
     * handles gradual underflow correctly.
     */
    if (biased <= 0) {
        return ml_ldexp_pure(mant_approx, (int)xi);
    }

    if (biased >= 0x7FF) {
        return ml_make_inf(0);
    }

    uint64_t bits = ((uint64_t)biased << 52);
    double res;
    memcpy(&res, &bits, sizeof(double));

    return res * mant_approx;
}

#endif /* LIBMATHC_FAST_MATH_H */
