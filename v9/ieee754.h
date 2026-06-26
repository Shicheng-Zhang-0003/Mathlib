#ifndef LIBMATHC_IEEE754_H
#define LIBMATHC_IEEE754_H

#include <stdint.h>
#include <math.h>

// Union for type-punning between double and 64-bit integer
typedef union { double d; uint64_t i; } ml_double_cast;

#define ML_LN2 0.693147180559945309417

// Pure IEEE 754 Bit-Masking Logarithm
inline double logarithm_ieee754(double x) {
    if (x <= 0.0) return 0.0 / 0.0;
    if (x == 1.0) return 0.0;

    ml_double_cast cast;
    cast.d = x;

    // Extract exponent (bits 52-62) and subtract bias (1023)
    int e = ((cast.i >> 52) & 0x7FF) - 1023;

    // Extract mantissa and restore the hidden bit (bit 52)
    uint64_t mantissa = (cast.i & 0xFFFFFFFFFFFFFULL) | 0x10000000000000ULL;

    // Convert mantissa to double in range [1.0, 2.0)
    double m = (double)mantissa / 4503599627370496.0; // 2^52

    // Adjust to [sqrt(2)/2, sqrt(2)] for optimal series convergence
    if (m > 1.4142135623730950) {
        m /= 2.0;
        e++;
    }

    // Fast series: 2 * (z + z^3/3 + z^5/5...) where z = (m-1)/(m+1)
    double z = (m - 1.0) / (m + 1.0);
    double z2 = z * z;
    double res = z;
    double term = z;

    for (int i = 3; i <= 15; i += 2) {
        term *= z2;
        res += term / i;
    }

    return 2.0 * res + e * ML_LN2;
}

// Pure IEEE 754 Bit-Masking Exponential
inline double exponential_ieee754(double x) {
    if (x == 0.0) return 1.0;

    // Calculate n = round(x / ln2)
    double n_d = x / ML_LN2;
    int n = (int)(n_d + (n_d > 0 ? 0.5 : -0.5));

    // Calculate remainder r = x - n * ln2
    double r = x - n * ML_LN2;

    // Taylor series for e^r (r is very small, converges instantly)
    double res = 1.0;
    double term = 1.0;
    for (int i = 1; i <= 15; i++) {
        term *= r / i;
        res += term;
    }

    // Multiply by 2^n using pure IEEE 754 exponent bit-shifting
    ml_double_cast cast;
    cast.d = res;
    cast.i += ((uint64_t)n << 52); // Add n to the exponent bits

    return cast.d;
}


// --- Pure Bitwise frexp and ldexp (No Standard Library) ---
inline double ml_ldexp_pure(double x, int exp) {
    ml_double_cast cast; cast.d = x;
    cast.i += ((uint64_t)exp << 52);
    return cast.d;
}

inline double ml_frexp_pure(double x, int *exp) {
    ml_double_cast cast; cast.d = x;
    // Extract exponent, adjust bias to 1022 for [0.5, 1.0) range
    *exp = ((cast.i >> 52) & 0x7FF) - 1022;
    // Mask out exponent, set it to 1022 (0x3FE)
    cast.i = (cast.i & 0x800FFFFFFFFFFFFFULL) | 0x3FE0000000000000ULL;
    return cast.d;
}

#endif
