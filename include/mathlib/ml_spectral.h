#ifndef MATHLIB_ML_SPECTRAL_H
#define MATHLIB_ML_SPECTRAL_H
#include "ml_compiler.h"
#include "ml_core.h"
#include "ml_types.h"
#include "ml_tensor.h"
#ifdef __cplusplus
extern "C" {
#endif
ML_API ml_status_t ml_cg_solve(ml_tensor_view_t A, const double *b, double *x, int max_iter, double tol);
ML_API ml_status_t ml_gmres_solve(ml_tensor_view_t A, const double *b, double *x, int restart, int max_iter, double tol);
ML_API ml_status_t ml_power_iter(ml_tensor_view_t A, double *lambda, double *vec, int max_iter, double tol);
ML_API ml_status_t ml_svd_jacobi(ml_tensor_view_t A, double *svals, ml_tensor_view_t U, ml_tensor_view_t Vt, int max_sweeps);
ML_API ml_status_t ml_qr_iter_eig(ml_tensor_view_t A, double *evals_re, double *evals_im, int max_iter, double tol);
#ifdef __cplusplus
}
#endif
#endif
