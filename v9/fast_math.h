#ifndef LIBMATHC_FAST_MATH_H
#define LIBMATHC_FAST_MATH_H
#include "bitwise_fp.h"

// Quake III Fast Inverse Sqrt for 64-bit doubles
inline double ml_fast_rsqrt(double number) {
    ml_fp_cast c; c.d = number;
    // Magic constant for 64-bit double precision
    c.u = 0x5fe6ec85e7de30daULL - (c.u >> 1);
    double y = c.d;
    // One iteration of Newton-Raphson for 64-bit precision
    return y * (1.5 - (number * 0.5 * y * y));
}

// Fast Log2 using the integer-float isomorphism
inline double ml_fast_log2(double x) {
    ml_fp_cast c; c.d = x;
    // Extract exponent, adjust bias (1023), and add mantissa approximation
    double exp = (double)((c.u >> 52) & 0x7FF) - 1023.0;
    double mant = (double)(c.u & ML_FP_MANT_MASK) / 4503599627370496.0;
    return exp + mant; // Linear approximation of log2(1+m)
}

// Fast Exp2 using the reverse isomorphism
inline double ml_fast_exp2(double x) {
    double exp_int = (double)(long long)x;
    double mant_frac = x - exp_int;
    ml_fp_cast c;
    // Construct the double directly from integer parts
    c.u = ((uint64_t)(exp_int + 1023) << 52) | (uint64_t)(mant_frac * 4503599627370496.0);
    return c.d;
}
#endif
