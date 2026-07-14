#ifndef LIBMATHC_COMBINATORICS_H
#define LIBMATHC_COMBINATORICS_H
#include <stdint.h>

static inline uint64_t factorial(int x) {
    if (x < 0 || x > 20) return 18446744073709551615ULL; // UINT64_MAX sentinel // Hard cap: 21! overflows uint64_t
    if (x == 0) return 1;
    uint64_t result = (uint64_t)x;
    while ((x - 1) > 0) {
        result *= (x - 1);
        x -= 1;
    }
    return result;
}

static inline uint64_t npr(int n, int r) {
    if (r < 0 || r > n) return 0;
    uint64_t num = factorial(n);
    uint64_t den = factorial(n - r);
    if (num == 18446744073709551615ULL || den == 18446744073709551615ULL || den == 0) return 18446744073709551615ULL;
    return num / den;
}

static inline uint64_t ncr(int n, int r) {
    if (r < 0 || r > n) return 0;
    uint64_t num = factorial(n);
    uint64_t den1 = factorial(n - r);
    uint64_t den2 = factorial(r);
    if (num == 18446744073709551615ULL || den1 == 18446744073709551615ULL || den2 == 18446744073709551615ULL || den1 == 0 || den2 == 0) return 18446744073709551615ULL;
    return num / (den1 * den2);
}
#endif
