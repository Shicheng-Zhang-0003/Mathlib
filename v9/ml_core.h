#ifndef ML_CORE_H
#define ML_CORE_H

#include <stdint.h>

typedef union { double d; uint64_t u; } ml_fp_cast;

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
    if (y == 0.0) return 0.0/0.0;
    double q = x / y;
    long long qi = (long long)q;
    double rem = x - qi * y;
    if (ml_fabs(rem) >= ml_fabs(y)) {
        rem = rem - (rem > 0 ? y : -y);
    }
    return rem;
}

// Pure C Round (Zero libm dependency)
static inline double ml_round(double x) {
    return (x >= 0.0) ? (double)(long long)(x + 0.5) : (double)(long long)(x - 0.5);
}
#endif
