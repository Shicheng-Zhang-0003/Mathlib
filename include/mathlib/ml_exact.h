#ifndef MATHLIB_ML_EXACT_H
#define MATHLIB_ML_EXACT_H

#include "ml_compiler.h"
#include "ml_types.h"
#include "ml_interval.h"
#include <limits.h>
#include <stddef.h>
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
/* Dynamically sized exact integers. Initialize before first use and clear
 * exactly once when finished. These owning structs must not be copied with C
 * assignment; use the copy functions. All output objects must already be
 * initialized. Capacity grows to available memory. */
typedef struct {
    uint32_t *limb; /* little-endian base 2^32 */
    size_t size;
    size_t capacity;
    int sign; /* 0 nonnegative, 1 negative; zero is always nonnegative */
} ml_mpz_t;
ML_API void ml_mpz_init (ml_mpz_t *a);
ML_API void ml_mpz_clear (ml_mpz_t *a);
ML_API int ml_mpz_is_valid (const ml_mpz_t *a);
ML_API int ml_mpz_is_zero (const ml_mpz_t *a);
ML_API int ml_mpz_sign (const ml_mpz_t *a); /* -1, 0, +1; invalid -> 0 */
ML_API size_t ml_mpz_bit_length (const ml_mpz_t *a);
ML_API int ml_mpz_compare (const ml_mpz_t *a, const ml_mpz_t *b); /* invalid -> INT_MIN */
ML_API ml_status_t ml_mpz_set_i64 (ml_mpz_t *r, int64_t value);
ML_API ml_status_t ml_mpz_set_u64 (ml_mpz_t *r, uint64_t value);
ML_API ml_status_t ml_mpz_set_string (ml_mpz_t *r, const char *text, unsigned base);
ML_API int ml_mpz_to_string (const ml_mpz_t *a, char *buffer, size_t size,
                             unsigned base); /* chars written; -1 on error */
ML_API ml_status_t ml_mpz_copy (ml_mpz_t *r, const ml_mpz_t *a);
ML_API ml_status_t ml_mpz_add (ml_mpz_t *r, const ml_mpz_t *a, const ml_mpz_t *b);
ML_API ml_status_t ml_mpz_sub (ml_mpz_t *r, const ml_mpz_t *a, const ml_mpz_t *b);
ML_API ml_status_t ml_mpz_mul (ml_mpz_t *r, const ml_mpz_t *a, const ml_mpz_t *b);
ML_API ml_status_t ml_mpz_divrem (ml_mpz_t *quotient, ml_mpz_t *remainder,
                                  const ml_mpz_t *a, const ml_mpz_t *b);
ML_API ml_status_t ml_mpz_gcd (ml_mpz_t *r, const ml_mpz_t *a, const ml_mpz_t *b);
ML_API ml_status_t ml_mpz_shift_left (ml_mpz_t *r, const ml_mpz_t *a, size_t bits);
ML_API double ml_mpz_to_double (const ml_mpz_t *a); /* rounded approximation */
/* Canonical arbitrary-precision rationals. The denominator is positive and
 * numerator/denominator are coprime. Initialize/clear and do not shallow-copy. */
typedef struct {
    ml_mpz_t numerator;
    ml_mpz_t denominator;
} ml_mpq_t;
ML_API void ml_mpq_init (ml_mpq_t *q);
ML_API void ml_mpq_clear (ml_mpq_t *q);
ML_API int ml_mpq_is_valid (const ml_mpq_t *q);
ML_API ml_status_t ml_mpq_set_i64 (ml_mpq_t *r, int64_t value);
ML_API ml_status_t ml_mpq_set_ratio (ml_mpq_t *r, const ml_mpz_t *numerator,
                                     const ml_mpz_t *denominator);
/* Parse a finite decimal literal exactly, e.g. "-12.50e-2". Accepts an
 * optional sign, decimal point, and base-10 exponent; no whitespace or
 * separators. At most 10000 significand digits and exponent magnitude 10000. */
ML_API ml_status_t ml_mpq_set_decimal (ml_mpq_t *r, const char *text);
ML_API ml_status_t ml_mpq_copy (ml_mpq_t *r, const ml_mpq_t *q);
ML_API ml_status_t ml_mpq_add (ml_mpq_t *r, const ml_mpq_t *a, const ml_mpq_t *b);
ML_API ml_status_t ml_mpq_sub (ml_mpq_t *r, const ml_mpq_t *a, const ml_mpq_t *b);
ML_API ml_status_t ml_mpq_mul (ml_mpq_t *r, const ml_mpq_t *a, const ml_mpq_t *b);
ML_API ml_status_t ml_mpq_div (ml_mpq_t *r, const ml_mpq_t *a, const ml_mpq_t *b);
ML_API ml_status_t ml_mpq_compare (const ml_mpq_t *a, const ml_mpq_t *b,
                                   int *comparison);
/* Produces the tightest adjacent-binary64 enclosure of the exact rational.
 * Returns ML_ERR_OVERFLOW when it lies outside the finite binary64 range. */
ML_API ml_status_t ml_mpq_to_interval (const ml_mpq_t *q, ml_interval_t *out);

#ifdef __cplusplus
}
#endif
#endif /* MATHLIB_ML_EXACT_H */
