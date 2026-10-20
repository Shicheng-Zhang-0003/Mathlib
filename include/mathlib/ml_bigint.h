#ifndef MATHLIB_ML_BIGINT_H
#define MATHLIB_ML_BIGINT_H

#include "ml_compiler.h"
#include "ml_core.h"
#include "ml_types.h"
#include <stdint.h>
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * v12R3: arbitrary-precision signed integers.
 *
 * Fixed-capacity, zero-allocation limb representation (base 2^32). Capacity
 * is 128 limbs = 4096 bits (~1233 decimal digits). Arithmetic that exceeds
 * this capacity returns ML_ERR_OVERFLOW. Values passed to operations must be
 * initialized with ml_bigint_init or another library operation; the public
 * representation is required to be canonical.
 *
 * Division truncates the quotient toward zero and gives the remainder the
 * dividend's sign. Modulo uses a non-negative Euclidean remainder with |b|.
 * ========================================================================== */

#define ML_BIGINT_LIMBS 128
typedef struct {
    uint32_t limb [ML_BIGINT_LIMBS];
    int n;    /* number of significant limbs; 0 means zero */
    int sign; /* 0 = >= 0, 1 = < 0 (zero is always sign 0) */
} ml_bigint_t;
ML_API void ml_bigint_init (ml_bigint_t *a);
ML_API ml_status_t ml_bigint_set_i64 (ml_bigint_t *a, int64_t v);
ML_API ml_status_t ml_bigint_set_u64 (ml_bigint_t *a, uint64_t v);
ML_API ml_status_t ml_bigint_copy (ml_bigint_t *r, const ml_bigint_t *a);
ML_API ml_status_t ml_bigint_from_str (ml_bigint_t *a, const char *s, int base);
ML_API int ml_bigint_to_str (const ml_bigint_t *a, char *buf, size_t bufsz, int base);

ML_API int ml_bigint_cmp (const ml_bigint_t *a, const ml_bigint_t *b);
ML_API int ml_bigint_is_zero (const ml_bigint_t *a);
ML_API int ml_bigint_bitlen (const ml_bigint_t *a);
ML_API double ml_bigint_to_double (const ml_bigint_t *a);

ML_API ml_status_t ml_bigint_add (ml_bigint_t *r, const ml_bigint_t *a, const ml_bigint_t *b);
ML_API ml_status_t ml_bigint_sub (ml_bigint_t *r, const ml_bigint_t *a, const ml_bigint_t *b);
ML_API ml_status_t ml_bigint_mul (ml_bigint_t *r, const ml_bigint_t *a, const ml_bigint_t *b);
ML_API ml_status_t ml_bigint_div (ml_bigint_t *q, ml_bigint_t *rem,
                                  const ml_bigint_t *a, const ml_bigint_t *b);
ML_API ml_status_t ml_bigint_mod (ml_bigint_t *r, const ml_bigint_t *a, const ml_bigint_t *b);
ML_API ml_status_t ml_bigint_pow (ml_bigint_t *r, const ml_bigint_t *base, uint64_t exp);
ML_API ml_status_t ml_bigint_powmod (ml_bigint_t *r, const ml_bigint_t *base,
                                     const ml_bigint_t *exp, const ml_bigint_t *mod);
ML_API ml_status_t ml_bigint_gcd (ml_bigint_t *r, const ml_bigint_t *a, const ml_bigint_t *b);
ML_API ml_status_t ml_bigint_lcm (ml_bigint_t *r, const ml_bigint_t *a, const ml_bigint_t *b);
ML_API ml_status_t ml_bigint_modinv (ml_bigint_t *r, const ml_bigint_t *a, const ml_bigint_t *m);
ML_API ml_status_t ml_bigint_shl (ml_bigint_t *r, const ml_bigint_t *a, unsigned bits);
ML_API ml_status_t ml_bigint_shr (ml_bigint_t *r, const ml_bigint_t *a, unsigned bits);
ML_API ml_status_t ml_bigint_abs (ml_bigint_t *r, const ml_bigint_t *a);
ML_API ml_status_t ml_bigint_neg (ml_bigint_t *r, const ml_bigint_t *a);

ML_API ml_status_t ml_bigint_factorial (ml_bigint_t *r, unsigned n);
ML_API ml_status_t ml_bigint_fib (ml_bigint_t *r, unsigned n);
ML_API ml_status_t ml_bigint_rand (ml_bigint_t *r, unsigned bits, uint64_t *state);
/* Deterministic Miller-Rabin probable-prime test. For inputs above 32 bits,
 * state must be non-null and rounds must be in [1,256]. Returns 1 for a
 * probable prime, 0 for composite or invalid input. The xorshift generator
 * used here is not cryptographically secure. */
ML_API int ml_bigint_is_probable_prime (const ml_bigint_t *a, int rounds,
                                        uint64_t *state);

#ifdef __cplusplus
}
#endif
#endif /* MATHLIB_ML_BIGINT_H */
