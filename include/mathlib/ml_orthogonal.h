#ifndef MATHLIB_ML_ORTHOGONAL_H
#define MATHLIB_ML_ORTHOGONAL_H

#include "ml_compiler.h"
#include "ml_core.h"

/* ============================================================================
 * Undergraduate orthogonal polynomials (Sturm-Liouville / special functions).
 *
 * Theorems covered: three-term recurrence, Rodrigues, orthogonality on
 * (-1,1) with weight 1 (Legendre), 1/sqrt(1-x^2) (Chebyshev T),
 * sqrt(1-x^2) (Chebyshev U), e^{-x^2} (Hermite phys), e^{-x} (Laguerre).
 *
 * Accuracy: all recurrences evaluated in long double (64-bit mantissa)
 * then rounded once to double. For n<=64, |x|<=2 the rounding error
 * dominates (<1 ULP vs mpmath on tested grid). Large n*x may grow to
 * ~n*eps; still backward-stable by three-term theory.
 * ========================================================================== */

ML_API double ml_legendre_p(int n, double x);
ML_API double ml_legendre_p_deriv(int n, double x);
ML_API double ml_chebyshev_t(int n, double x);
ML_API double ml_chebyshev_u(int n, double x);
ML_API double ml_hermite_h(int n, double x);
ML_API double ml_laguerre_l(int n, double x);
ML_API double ml_laguerre_l_gen(int n, double alpha, double x);

#endif /* MATHLIB_ML_ORTHOGONAL_H */
