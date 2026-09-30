#ifndef MATHLIB_ML_POLYNOMIAL_H
#define MATHLIB_ML_POLYNOMIAL_H

#include "ml_compiler.h"
#include "ml_core.h"
#ifdef __cplusplus
extern "C" {
#endif

/* FIX: Added const to non-mutating coefficient arrays */
/* Polynomial contract:
 * - coeffs[0..degree] with degree >= 0; caller provides degree+1 entries.
 *   No explicit degree cap is enforced, but very large degree may overflow
 *   to Inf and Horner cost is O(degree).
 * - No deflation: ml_polynomial_newton finds a single root near x0; it does
 *   not deflate or enumerate all roots.
 * - Finite-input validation: non-finite coeffs/x/x0/epsilon -> NaN (eval and
 *   newton return NaN; derivative fills outputs with NaN).
 * - Convergence-failure signal: newton returns NaN when max_iter is
 *   exhausted, the derivative is exactly zero, or an iterate is non-finite.
 */
ML_API double ml_polynomial_eval(const double *coeffs, int degree, double x);
ML_API void ml_polynomial_derivative(const double *coeffs, int degree, double *out);
ML_API double ml_polynomial_newton(const double *coeffs, int degree, double x0, double epsilon, int max_iter);

#ifdef __cplusplus
}
#endif
#endif /* MATHLIB_ML_POLYNOMIAL_H */
