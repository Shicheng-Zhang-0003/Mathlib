#ifndef LIBMATHC_COMBINATORICS_H
#define LIBMATHC_COMBINATORICS_H

//Library header file for combinatorics
#include <stdio.h>
static inline uint64_t factorial(int x);
static inline uint64_t factorial(int x) {
    if (x < 0) {return (uint64_t) (0);}
    if (x == 0) {return (uint64_t) (1);}
    uint64_t result = (uint64_t) (x);
    while ((x - 1) > 0) {result *= (x - 1);
    x -= 1;} return result;
} static inline uint64_t npr(int n, int r);
static inline uint64_t npr(int n, int r) {
    if (r < 0 || r > n) {return (uint64_t) (0);}
    return factorial (n) / factorial (n - r);
}
static inline uint64_t ncr(int n, int r);
static inline uint64_t ncr(int n, int r) {
    if (r < 0 || r > n) {return (uint64_t) (0);}
    return factorial (n) / (factorial (n - r) * factorial (r));
}



#endif
