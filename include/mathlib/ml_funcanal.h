#ifndef MATHLIB_ML_FUNCANAL_H
#define MATHLIB_ML_FUNCANAL_H

#include "ml_compiler.h"
#include "ml_core.h"
#include "ml_types.h"
#include "ml_numerical.h"
#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Functional analysis, computable slice: Banach contraction solver with a
 * caller-certified Lipschitz constant and an a posteriori bound, Neumann
 * series inverses, Hilbert-space best approximation (MGS + back
 * substitution), and dominant-eigenvalue Rayleigh quotients (wrapping the
 * staged power iteration).
 *
 * No oracle ULP claimed (regression-guard tier). Caps keep frames < 8KB.
 * ========================================================================== */

/* Banach: x_{n+1} = T(x_n); L in [0,1) is the caller certificate.
 * bound_out (optional) receives L^n/(1-L)|x_1-x_0|. */
ML_API ml_status_t ml_banach_fixed_point(ml_func_t T, double x0, double L,
                                         double tol, int max_iter,
                                         double *xstar, double *bound_out);
/* Neumann: out = sum_{k<terms} (I-A)^k, needs row-sum ||I-A|| < 1, n <= 16. */
ML_API ml_status_t ml_neumann_inverse(const double *A, int n, int terms, double *out);
/* Best approximation: coeff minimizing ||V coeff - y||, V m x k row-major
 * (m <= 1024, k <= 32); SINGULAR on rank loss. */
ML_API ml_status_t ml_hilbert_project(const double *V, const double *y,
                                      int m, int k, double *coeff);
/* Largest-magnitude eigenvalue of symmetric A (row-major, n <= 128). */
ML_API ml_status_t ml_rayleigh_max(const double *A, int n, double *lambda);

#ifdef __cplusplus
}
#endif
#endif /* MATHLIB_ML_FUNCANAL_H */
