#ifndef MATHLIB_ML_POLYQ_H
#define MATHLIB_ML_POLYQ_H

#include "ml_compiler.h"
#include "ml_exact.h"
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
/* Exact univariate polynomials over Q. Coefficients are stored in ascending
 * degree order. The owning struct must be initialized/cleared and must not be
 * shallow-copied. Capacity and running time are limited by available memory. */
typedef struct {
    ml_mpq_t *coefficient;
    size_t size; /* active coefficient count; zero polynomial has size 0 */
    size_t capacity;
} ml_polyq_t;
ML_API void ml_polyq_init (ml_polyq_t *p);
ML_API void ml_polyq_clear (ml_polyq_t *p);
ML_API int ml_polyq_is_valid (const ml_polyq_t *p);
ML_API int ml_polyq_is_zero (const ml_polyq_t *p);
ML_API ptrdiff_t ml_polyq_degree (const ml_polyq_t *p); /* zero polynomial: -1 */
ML_API ml_status_t ml_polyq_copy (ml_polyq_t *r, const ml_polyq_t *p);
ML_API ml_status_t ml_polyq_set_coefficient (ml_polyq_t *p, size_t degree,
                                             const ml_mpq_t *value);
ML_API ml_status_t ml_polyq_get_coefficient (const ml_polyq_t *p, size_t degree,
                                             ml_mpq_t *value);
ML_API ml_status_t ml_polyq_equal (const ml_polyq_t *a, const ml_polyq_t *b,
                                   int *equal);
ML_API ml_status_t ml_polyq_add (ml_polyq_t *r, const ml_polyq_t *a,
                                 const ml_polyq_t *b);
ML_API ml_status_t ml_polyq_sub (ml_polyq_t *r, const ml_polyq_t *a,
                                 const ml_polyq_t *b);
ML_API ml_status_t ml_polyq_mul (ml_polyq_t *r, const ml_polyq_t *a,
                                 const ml_polyq_t *b);
ML_API ml_status_t ml_polyq_derivative (ml_polyq_t *r, const ml_polyq_t *p);
ML_API ml_status_t ml_polyq_evaluate (ml_mpq_t *value, const ml_polyq_t *p,
                                      const ml_mpq_t *x);
ML_API ml_status_t ml_polyq_divrem (ml_polyq_t *quotient,
                                    ml_polyq_t *remainder,
                                    const ml_polyq_t *dividend,
                                    const ml_polyq_t *divisor);
ML_API ml_status_t ml_polyq_gcd (ml_polyq_t *gcd, const ml_polyq_t *a,
                                 const ml_polyq_t *b);
/* Verifies the complete exact witness: divisor is nonzero,
 * dividend = divisor*quotient+remainder, and deg(remainder)<deg(divisor). */
ML_API ml_status_t ml_polyq_verify_division (const ml_polyq_t *dividend,
                                             const ml_polyq_t *divisor,
                                             const ml_polyq_t *quotient,
                                             const ml_polyq_t *remainder,
                                             int *valid);

#ifdef __cplusplus
}
#endif
#endif /* MATHLIB_ML_POLYQ_H */
