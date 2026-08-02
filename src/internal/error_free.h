#ifndef LIBMATHC_ERROR_FREE_H
#define LIBMATHC_ERROR_FREE_H

/* MATHLIB_V12A1_ERROR_FREE_CLEANUP */
/*
 * ERROR-FREE TRANSFORMATIONS
 *
 * These functions compute exact rounding errors from floating-point
 * operations. They are building blocks for compensated algorithms
 * (Kahan summation, compensated Horner, etc.).
 *
 * NAMING CONVENTION:
 *   ml_two_sum(a, b, &err)      -> s = fl(a + b), err = (a + b) - s exactly
 *   ml_fast_two_sum(a, b, &err) -> same, but REQUIRES |a| >= |b|
 *   ml_two_product(a, b, &err)  -> p = fl(a * b), err = (a * b) - p exactly
 *   ml_fma_soft(a, b, c)        -> software FMA via Two-Product + Two-Sum
 *
 * IMPORTANT: ml_fma_soft is NOT the same as ML_FMA (from ml_compiler.h).
 *
 *   ML_FMA(a, b, c)    = hardware fused multiply-add, rounds ONCE (0.5 ULP)
 *   ml_fma_soft(a, b, c) = software emulation, rounds TWICE (up to 2 ULP)
 *
 * Use ML_FMA for all production math. Use ml_fma_soft only when you
 * explicitly need the software path (e.g., testing, or platforms
 * without hardware FMA).
 */

#include "ml_compiler.h"

/*
 * Knuth's Two-Sum.
 * No magnitude assumption on a, b.
 * Returns s = fl(a + b) and sets *err = (a + b) - s exactly.
 */
static inline double ml_two_sum(double a, double b, double *err) {
    double s = a + b;
    double v = s - a;
    *err = (a - (s - v)) + (b - v);
    return s;
}

/*
 * Dekker's Fast Two-Sum.
 * PRECONDITION: |a| >= |b| (caller must guarantee this).
 * Returns s = fl(a + b) and sets *err = (a + b) - s exactly.
 *
 * Faster than ml_two_sum (3 ops vs 6) but only valid when
 * the magnitude ordering is known.
 */
static inline double ml_fast_two_sum(double a, double b, double *err) {
    double s = a + b;
    double z = s - a;
    *err = b - z;
    return s;
}

/*
 * Dekker's Two-Product.
 * Returns p = fl(a * b) and sets *err = (a * b) - p exactly.
 *
 * Uses Dekker splitting (multiply by 2^26 + 1) to split each
 * operand into high and low 26-bit halves. This is exact for
 * any finite double because the significand is 53 bits.
 *
 * Note: If hardware FMA is available, you can compute the error
 * term more cheaply as: err = ML_FMA(a, b, -p). But this function
 * remains useful for platforms without FMA and for code that
 * needs to be FMA-independent.
 */
static inline double ml_two_product(double a, double b, double *err) {
    double p = a * b;
    double ca = a * 67108865.0; /* 2^26 + 1 */
    double a_hi = ca - (ca - a);
    double a_lo = a - a_hi;
    double cb = b * 67108865.0;
    double b_hi = cb - (cb - b);
    double b_lo = b - b_hi;
    *err = ((a_hi * b_hi - p) + a_hi * b_lo + a_lo * b_hi) + a_lo * b_lo;
    return p;
}

/*
 * Software FMA emulation.
 *
 * Computes fl(a * b + c) using error-free transformations.
 * This rounds TWICE and has up to 2 ULP error.
 *
 * DO NOT USE THIS IN PRODUCTION MATH. Use ML_FMA instead.
 * This exists for:
 *   - testing error-free transformation correctness
 *   - platforms without hardware FMA (the ML_FMA macro falls
 *     back to (a*b)+c on such platforms, which is even worse)
 *   - educational/reference purposes
 */
static inline double ml_fma_soft(double a, double b, double c) {
    double p, prod_err;
    p = ml_two_product(a, b, &prod_err);
    double s1, sum_err;
    s1 = ml_two_sum(p, c, &sum_err);
    return s1 + (prod_err + sum_err);
}

#endif /* LIBMATHC_ERROR_FREE_H */
