#ifndef MATHLIB_ML_NUMBERTHEORY_H
#define MATHLIB_ML_NUMBERTHEORY_H

#include "ml_compiler.h"
#include "ml_core.h"
#include "ml_types.h"
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Undergraduate number theory: divisibility, congruences, reciprocity.
 *
 * Theorems: Euler criterion, quadratic reciprocity, Chinese remainder
 * (generalized), Fermat/Euler/Miller-Rabin, Mobius inversion, Dirichlet
 * convolution identities, primitive-root structure of (Z/pZ)*.
 *
 * Accuracy: exact integer arithmetic (0 ULP); returns UINT64_MAX / 2
 * sentinels only when input exceeds 64-bit exact range.
 * Mobius contract: -1/0/1 for decided inputs, 2 = UNRESOLVED when the
 * cofactor beyond the 1e6 trial bound is composite with unknown
 * square-freeness. Treat 2 as unknown, never as a Mobius value.
 * crt2 contract: generalized 2-congruence; returns SINGULAR when
 * incompatible (d1 % gcd != 0), WORKSPACE on u64 overflow of the lcm.
 * mult_order contract: exact order (divisor of phi); the factorizer
 * strips small factors by trial then reduces any residual cofactor via
 * gcd-guarded division, so a composite residual cannot over-divide.
 * ========================================================================== */

ML_API int ml_legendre_symbol (int64_t a, uint64_t p);
ML_API int ml_jacobi_symbol (int64_t a, uint64_t n);
ML_API int ml_kronecker_symbol (int64_t a, int64_t n);
ML_API int ml_mobius (int64_t n);
ML_API uint64_t ml_mult_order (uint64_t a, uint64_t m);
ML_API int ml_is_primitive_root (uint64_t a, uint64_t p);
ML_API ml_status_t ml_crt2 (uint64_t a1, uint64_t m1, uint64_t a2, uint64_t m2,
                            uint64_t *x, uint64_t *lcm);
ML_API uint64_t ml_prime_pi (uint64_t n);
ML_API int ml_quadratic_reciprocity_check (uint64_t p, uint64_t q);

/* ---- Olympiad kit (exact u64, UINT64_MAX sentinel) ---- */
/* Legendre: v_p(n!) for prime p (non-prime p -> UINT64_MAX). */
ML_API uint64_t ml_vp_factorial (long long n, uint64_t p);
/* AtCoder floor_sum: sum_{0<=i<n} floor((a*i+b)/m); m = 0 -> UINT64_MAX. */
ML_API uint64_t ml_floor_sum_u64 (uint64_t n, uint64_t m, uint64_t a, uint64_t b);
/* LTE: v_p(a^n-b^n). Odd p needs p | a-b, p coprime ab; p = 2 needs a,b
 * odd (n odd -> v2(a-b); n even -> v2(a-b)+v2(a+b)+v2(n)-1). */
ML_API uint64_t ml_lte_vp (uint64_t a, uint64_t b, uint64_t n, uint64_t p);
/* sqrt(D) continued-fraction period: >=1 period, 0 = perfect square,
 * -1 = invalid (D<2), -2 = period exceeds the 10000 bound. */
ML_API int ml_sqrt_cf_period (uint64_t D);
/* k-th convergent (0-indexed) of sqrt(D); OVERFLOW when u64 overflows. */
ML_API ml_status_t ml_cf_convergent_sqrt (uint64_t D, int k, uint64_t *num, uint64_t *den);
/* Minimal x^2-Dy^2 = 1. Square D: INVALID_ARG; solution beyond u64:
 * OVERFLOW; period beyond 20000 convergents: UNSUPPORTED. */
ML_API ml_status_t ml_pell_fund (uint64_t D, uint64_t *x, uint64_t *y);

#ifdef __cplusplus
}
#endif
#endif /* MATHLIB_ML_NUMBERTHEORY_H */
