#ifndef LIBMATHC_FAST_MATH_H
#define LIBMATHC_FAST_MATH_H
#include "bitwise_fp.h"

// Quake III Fast Inverse Sqrt for 64-bit doubles
static inline double ml_fast_rsqrt(double number) {
    ml_fp_cast c; c.d = number;
    // Magic constant for 64-bit double precision
    c.u = 0x5fe6ec85e7de30daULL - (c.u >> 1);
    double y = c.d;
    // One iteration of Newton-Raphson for 64-bit precision
    return y * (1.5 - (number * 0.5 * y * y));
}

// Fast Log2 using the integer-float isomorphism
static inline double ml_fast_log2(double x) {
    if (x < 0.0) return 0.0/0.0;
    if (x == 0.0) return -1.0/0.0;
    if (ml_isinf(x) || ml_isnan(x)) return x;
    ml_fp_cast c; c.d = x;
    // Extract exponent, adjust bias (1023), and add mantissa approximation
    double exp = (double)((c.u >> 52) & 0x7FF) - 1023.0;
    double mant = (double)(c.u & ML_FP_MANT_MASK) / 4503599627370496.0;
    return exp + mant; // Linear approximation of log2(1+m)
}

// Fast Exp2 using reverse isomorphism + Minimax fractional polynomial
static inline double ml_fast_exp2(double x) {
    long long xi = (long long)x;
    if (x < 0.0 && x != (double)xi) xi -= 1; // Bitwise floor for negative numbers
    double exp_int = (double)xi;
    double frac = x - exp_int;
    // Minimax polynomial for 2^frac on [0, 1]
    double p = frac * (0.6931471805599453 + frac * (0.2402265069591007 + frac * 0.05550410866482158));
    double mant_approx = 1.0 + p;
    ml_fp_cast c;
    c.u = ((uint64_t)((long long)exp_int + 1023) << 52);
    return c.d * mant_approx;
}
#endif
