#ifndef MATHLIB_ML_INTERVAL_H
#define MATHLIB_ML_INTERVAL_H

#include "ml_compiler.h"
#include "ml_types.h"
#ifdef __cplusplus
extern "C" {
#endif
/*
 * Outward-rounded real intervals for basic arithmetic.
 *
 * This layer deliberately does not wrap transcendental libm functions.
 * Certified arithmetic is supported only for IEEE-754 binary64 with the
 * runtime rounding mode set to round-to-nearest. Each arithmetic endpoint is
 * widened by one representable value, so the exact real result is enclosed
 * under those platform/compiler assumptions. Inputs and returned endpoints
 * must be finite; an operation whose enclosure would overflow returns
 * ML_ERR_OVERFLOW and leaves its output unchanged.
 */
typedef struct {
    double lo;
    double hi;
} ml_interval_t;
ML_API int ml_interval_is_supported (void);
ML_API int ml_interval_is_valid (ml_interval_t x);
ML_API ml_status_t ml_interval_make (double lo, double hi, ml_interval_t *out);
ML_API ml_status_t ml_interval_point (double x, ml_interval_t *out);

ML_API ml_status_t ml_interval_add (ml_interval_t a, ml_interval_t b,
                                    ml_interval_t *out);
ML_API ml_status_t ml_interval_sub (ml_interval_t a, ml_interval_t b,
                                    ml_interval_t *out);
ML_API ml_status_t ml_interval_mul (ml_interval_t a, ml_interval_t b,
                                    ml_interval_t *out);
/* Division is undefined for an interval denominator containing zero. */
ML_API ml_status_t ml_interval_div (ml_interval_t a, ml_interval_t b,
                                    ml_interval_t *out);
ML_API ml_status_t ml_interval_neg (ml_interval_t a, ml_interval_t *out);
ML_API ml_status_t ml_interval_hull (ml_interval_t a, ml_interval_t b,
                                     ml_interval_t *out);
ML_API int ml_interval_contains (ml_interval_t x, double value);
ML_API int ml_interval_contains_zero (ml_interval_t x);

#ifdef __cplusplus
}
#endif
#endif /* MATHLIB_ML_INTERVAL_H */
