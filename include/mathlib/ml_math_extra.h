#ifndef MATHLIB_ML_MATH_EXTRA_H
#define MATHLIB_ML_MATH_EXTRA_H

#include "ml_compiler.h"
#include "ml_core.h"
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * v12R3: libm ABI parity layer.
 *
 * Supplementary functions corresponding to operations in ISO C <math.h>.
 * Their specific edge behavior is documented below; this layer does not
 * claim universal correctly-rounded equivalence with a platform libm.
 * ========================================================================== */

/* Library-local classification codes returned by ml_fpclassify. */
#define ML_CLASS_NAN 0
#define ML_CLASS_INFINITE 1
#define ML_CLASS_NORMAL 2
#define ML_CLASS_SUBNORMAL 3
#define ML_CLASS_ZERO 4

/* Adjacent floating-point navigation */
ML_API double ml_nextafter (double x, double y);
ML_API double ml_nexttoward (double x, long double y);

/* Significand / exponent manipulation */
ML_API double ml_modf (double x, double *iptr);
ML_API double ml_scalbn (double x, int n);
ML_API double ml_scalbln (double x, long n);
ML_API int ml_ilogb (double x);
ML_API double ml_logb_fp (double x);

/* Rounding (nearest-even) and conversions. Integer conversions return the
 * matching MIN value for NaN, infinity, or an out-of-range rounded value. */
ML_API double ml_rint (double x);
ML_API double ml_nearbyint (double x);
ML_API long ml_lrint (double x);
ML_API long long ml_llrint (double x);
ML_API long ml_lround (double x);
ML_API long long ml_llround (double x);

/* Remainder with signed low seven quotient bits, classification, FMA. */
ML_API double ml_remquo (double x, double y, int *quo);
ML_API int ml_fpclassify (double x);
ML_API int ml_isnormal (double x);
ML_API double ml_fma (double a, double b, double c);

/* Convenience: simultaneous sine and cosine, specials, tgamma sign */
ML_API void ml_sincos (double x, double *s, double *c);
ML_API double ml_nan (const char *tag);
ML_API double ml_inf (void);
ML_API double ml_tgamma (double x);
/* Returns log(|Gamma(x)|). signgamp is -1/+1 for finite non-poles, 0 for
 * NaN, negative infinity, and finite negative-integer poles. */
ML_API double ml_lgamma_r (double x, int *signgamp);

#ifdef __cplusplus
}
#endif
#endif /* MATHLIB_ML_MATH_EXTRA_H */
