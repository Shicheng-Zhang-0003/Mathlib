#ifndef MATHLIB_V10_LINALG_H
#define MATHLIB_V10_LINALG_H

#include "ml_compiler.h"
#include "ml_tensor.h"
#include "ml_core.h"
#include "ml_types.h"

ML_API ml_status_t ml_lu_decomp(ml_tensor_view_t A, ml_tensor_view_t LU, int* P, ml_workspace_t* ws);
ML_API ml_status_t ml_solve(ml_tensor_view_t A, double* b, double* x, ml_workspace_t* ws);
ML_API ml_status_t ml_cholesky(ml_tensor_view_t A, ml_tensor_view_t L);
ML_API ml_status_t ml_qr_solve(ml_tensor_view_t A, const double* b, double* x, ml_workspace_t* ws);
ML_API ml_status_t ml_solve_refined(ml_tensor_view_t A, double* b, double* x, ml_workspace_t* ws);
ML_API double ml_determinant(ml_tensor_view_t A, ml_workspace_t* ws);
ML_API ml_status_t ml_inverse(ml_tensor_view_t A, ml_tensor_view_t Inv, ml_workspace_t* ws);
ML_API void ml_transpose(ml_tensor_view_t A, ml_tensor_view_t T);
ML_API double ml_dot(const double* a, const double* b, int n);
ML_API void ml_cross3(const double a[3], const double b[3], double out[3]);
ML_API double ml_norm(const double* a, int n);
ML_API double ml_point_line_dist2d(double px, double py, double ax, double ay, double bx, double by);
ML_API double ml_point_plane_dist(double px, double py, double pz, double a, double b, double c, double d);
ML_API ml_status_t ml_eigen2x2(double a, double b, double c, double d, double *l0, double *l1);
/* Undergraduate spectral theorems: Jacobi symmetric eigensolver
 * (spectral theorem), analytic SVD 2x2 (SVD existence), 2x2 matrix
 * exponential (functional calculus / Jordan form). Jacobi is backward
 * stable; 2x2 closed forms are <1 ULP. */
ML_API ml_status_t ml_jacobi_eigen_symmetric(ml_tensor_view_t A, double *evals,
                                             ml_tensor_view_t V, int max_sweeps);
ML_API ml_status_t ml_svd_2x2(double a, double b, double c, double d,
                              double *s0, double *s1);
ML_API ml_status_t ml_matrix_exp_2x2(double a, double b, double c, double d,
                                     double *e00, double *e01,
                                     double *e10, double *e11);


/* Matrix-Vector Multiplication (y = Ax) */
ML_API void ml_matvec(ml_tensor_view_t A, const double* x, double* out);

#endif /* MATHLIB_V10_LINALG_H */
