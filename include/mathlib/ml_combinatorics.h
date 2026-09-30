#ifndef MATHLIB_ML_COMBINATORICS_H
#define MATHLIB_ML_COMBINATORICS_H

#include <stdint.h>
#include "ml_compiler.h"
#ifdef __cplusplus
extern "C" {
#endif

ML_API uint64_t ml_factorial(int x);
ML_API uint64_t ml_npr(int n, int r);
ML_API uint64_t ml_ncr(int n, int r);
ML_API uint64_t ml_gcd(uint64_t a, uint64_t b);
ML_API uint64_t ml_lcm(uint64_t a, uint64_t b);
ML_API uint64_t ml_modpow(uint64_t a, uint64_t e, uint64_t m);
ML_API int64_t ml_egcd(int64_t a, int64_t b, int64_t *x, int64_t *y);
ML_API uint64_t ml_modinv(uint64_t a, uint64_t m);
ML_API int ml_is_prime64(uint64_t n);
ML_API uint64_t ml_phi(uint64_t n);
ML_API uint64_t ml_catalan(int n);
ML_API uint64_t ml_fib_pair(long long n, uint64_t *fn, uint64_t *fn1);
ML_API uint64_t ml_derange(int n);

#ifdef __cplusplus
}
#endif
#endif /* MATHLIB_ML_COMBINATORICS_H */
