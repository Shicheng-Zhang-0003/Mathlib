#ifndef LIBMATHC_COMBINATORICS_H
#define LIBMATHC_COMBINATORICS_H
#include <stdint.h>

static inline uint64_t factorial(int x) {
    if (x < 0) return 0;
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
    return factorial(n) / factorial(n - r);
}

static inline uint64_t ncr(int n, int r) {
    if (r < 0 || r > n) return 0;
    return factorial(n) / (factorial(n - r) * factorial(r));
}
#endif
