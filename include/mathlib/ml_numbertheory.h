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

ML_API int ml_legendre_symbol(int64_t a, uint64_t p);
ML_API int ml_jacobi_symbol(int64_t a, uint64_t n);
ML_API int ml_kronecker_symbol(int64_t a, int64_t n);
ML_API int ml_mobius(int64_t n);
ML_API uint64_t ml_mult_order(uint64_t a, uint64_t m);
ML_API int ml_is_primitive_root(uint64_t a, uint64_t p);
ML_API ml_status_t ml_crt2(uint64_t a1, uint64_t m1, uint64_t a2, uint64_t m2,
                           uint64_t *x, uint64_t *lcm);
ML_API uint64_t ml_prime_pi(uint64_t n);
ML_API int ml_quadratic_reciprocity_check(uint64_t p, uint64_t q);

#ifdef __cplusplus
}
#endif
#endif /* MATHLIB_ML_NUMBERTHEORY_H */
