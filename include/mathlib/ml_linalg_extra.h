#ifndef MATHLIB_ML_LINALG_EXTRA_H
#define MATHLIB_ML_LINALG_EXTRA_H

#include "ml_compiler.h"
#include "ml_tensor.h"
#include "ml_core.h"
#include "ml_types.h"
#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * v12R2 LINEAR-ALGEBRA DEPTH PACK
 *
 * Graduate / research-grade decompositions, iterative solvers and matrix
 * functions that complement the undergraduate core in ml_linalg.h.
 *
 * Storage is row-major (see ml_tensor_view_t); all routines are stateless and
 * thread-safe. Routines that allocate scratch take an ml_workspace_t and
 * return ML_ERR_WORKSPACE when the caller's buffer is too small.
 * ========================================================================== */

/* Householder reduction of a general square matrix to upper Hessenberg form:
 * A = Q H Q^T. Q (n x n) may be NULL if the transform is not wanted. */
ML_API ml_status_t ml_hessenberg (ml_tensor_view_t A, ml_tensor_view_t H,
                                  ml_tensor_view_t Q);

/* Thin singular value decomposition A = U S V^T for a general m x n matrix.
 * S has k = min(m,n) entries sorted descending. U is m x k (pass .data==NULL
 * to skip), Vt is k x n (V^T; pass .data==NULL to skip). One-sided Jacobi;
 * high relative accuracy on well-scaled inputs. */
ML_API ml_status_t ml_svd (ml_tensor_view_t A, ml_tensor_view_t U, double *S,
                           ml_tensor_view_t Vt, ml_workspace_t *ws);

/* Moore-Penrose pseudoinverse computed from the SVD. Apinv is n x m. */
ML_API ml_status_t ml_pinv (ml_tensor_view_t A, ml_tensor_view_t Apinv,
                            double tol, ml_workspace_t *ws);

/* 2-norm condition number sigma_max / sigma_min (>= 1). HI for singular. */
ML_API ml_status_t ml_cond2 (ml_tensor_view_t A, double *out_cond,
                             ml_workspace_t *ws);

/* Numerical rank: number of singular values greater than tol. */
ML_API ml_status_t ml_matrix_rank (ml_tensor_view_t A, double tol, int *rank,
                                   ml_workspace_t *ws);

/* Iterative linear solvers. A is n x n; b, x are length n.
 * Return ML_SUCCESS on convergence, ML_ERR_CONVERGENCE otherwise.
 * *iters (optional) receives the iteration count. */
ML_API ml_status_t ml_gauss_seidel (ml_tensor_view_t A, const double *b,
                                    double *x, int max_iter, double tol,
                                    int *iters);
ML_API ml_status_t ml_sor (ml_tensor_view_t A, const double *b, double *x,
                           double omega, int max_iter, double tol, int *iters);
ML_API ml_status_t ml_conjugate_gradient (ml_tensor_view_t A, const double *b,
                                          double *x, int max_iter, double tol,
                                          int *iters);

/* Matrix functions (functional calculus) for a general square matrix.
 * ml_matrix_exp uses scaling-and-squaring with a Taylor kernel;
 * ml_matrix_sqrt uses the Denman-Beavers iteration; ml_matrix_log uses
 * inverse scaling-and-squaring; ml_matrix_pow(A,p) = exp(p log A). */
ML_API ml_status_t ml_matrix_exp (ml_tensor_view_t A, ml_tensor_view_t E,
                                  ml_workspace_t *ws);
ML_API ml_status_t ml_matrix_sqrt (ml_tensor_view_t A, ml_tensor_view_t R,
                                   ml_workspace_t *ws);
ML_API ml_status_t ml_matrix_log (ml_tensor_view_t A, ml_tensor_view_t L,
                                  ml_workspace_t *ws);
ML_API ml_status_t ml_matrix_pow (ml_tensor_view_t A, double p,
                                  ml_tensor_view_t P, ml_workspace_t *ws);

#ifdef __cplusplus
}
#endif
#endif /* MATHLIB_ML_LINALG_EXTRA_H */
