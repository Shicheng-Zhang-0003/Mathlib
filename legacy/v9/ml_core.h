#ifndef ML_CORE_H
#define ML_CORE_H

#include <stdint.h>

typedef union { double d; uint64_t u; } ml_fp_cast; // Union type-punning (Universally supported by GCC/Clang/MSVC)

// Pure Bitmask Absolute Value
static inline double ml_fabs(double x) {
    ml_fp_cast c; c.d = x;
    c.u &= 0x7FFFFFFFFFFFFFFFULL;
    return c.d;
}

// Pure Bitmask NaN Check
static inline int ml_isnan(double x) {
    ml_fp_cast c; c.d = x;
    uint64_t exp = c.u & 0x7FF0000000000000ULL;
    uint64_t mant = c.u & 0x000FFFFFFFFFFFFFULL;
    return (exp == 0x7FF0000000000000ULL) && (mant != 0);
}

// Pure Bitmask Infinity Check
static inline int ml_isinf(double x) {
    ml_fp_cast c; c.d = x;
    uint64_t exp = c.u & 0x7FF0000000000000ULL;
    uint64_t mant = c.u & 0x000FFFFFFFFFFFFFULL;
    return (exp == 0x7FF0000000000000ULL) && (mant == 0);
}

// Hardware SQRT via Inline Assembly (Zero libm dependency)
static inline double ml_sqrt(double x) {
    if (x < 0.0) return 0.0/0.0;
    if (x == 0.0) return 0.0;
    double res;
    __asm__ ("sqrtsd %1, %0" : "=x" (res) : "x" (x));
    return res;
}

// Pure C Software Fmod (Zero libm dependency)
static inline double ml_fmod(double x, double y) {
    if (ml_isnan(x) || ml_isnan(y) || ml_isinf(x)) return 0.0/0.0;
    if (ml_isinf(y)) return x;
    if (y == 0.0) return 0.0/0.0;
    double q = x / y;
    if (q > 9.22e18) q = 9.22e18;
    if (q < -9.22e18) q = -9.22e18;
    long long qi = (long long)q;
    double rem = x - qi * y;
    if (ml_fabs(rem) >= ml_fabs(y)) {
        rem = rem - (rem > 0 ? y : -y);
    }
    return rem;
}

// Pure C Round (Zero libm dependency)
static inline double ml_round(double x) {
    if (ml_isnan(x) || ml_isinf(x)) return x;
    if (x > 9.22e18) return x;
    if (x < -9.22e18) return x;
    return (x >= 0.0) ? (double)(long long)(x + 0.5) : (double)(long long)(x - 0.5);
}

// Pure Bitmask Copysign
static inline double ml_copysign(double x, double y) {
    ml_fp_cast cx, cy;
    cx.d = x; cy.d = y;
    cx.u = (cx.u & 0x7FFFFFFFFFFFFFFFULL) | (cy.u & 0x8000000000000000ULL);
    return cx.d;
}

// --- Pure Bitwise frexp and ldexp (No Standard Library) ---
static inline double ml_ldexp_pure(double x, int exp) {
    ml_fp_cast cast; cast.d = x;
    cast.u += ((uint64_t)exp << 52);
    return cast.d;
}

static inline double ml_frexp_pure(double x, int *exp) {
    ml_fp_cast cast; cast.d = x;
    // Extract exponent, adjust bias to 1022 for [0.5, 1.0) range
    *exp = ((cast.u >> 52) & 0x7FF) - 1022;
    // Mask out exponent, set it to 1022 (0x3FE)
    cast.u = (cast.u & 0x800FFFFFFFFFFFFFULL) | 0x3FE0000000000000ULL;
    return cast.d;
}


#endif
