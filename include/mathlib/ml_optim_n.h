#ifndef MATHLIB_ML_OPTIM_N_H
#define MATHLIB_ML_OPTIM_N_H
#include "ml_compiler.h"
#include "ml_core.h"
#include "ml_types.h"
typedef double (*ml_vec_func_t)(const double *x, int n);
/* Dimension caps (fixed stack workspaces, no heap): Nelder-Mead n<=16,
 * L-BFGS / Adam n<=32. Larger n returns ML_ERR_INVALID_ARG. */
/* ml_nelder_mead: converges when BOTH the function range |Fhi-Flo| and the
 * simplex diameter max||S_i-S_lo|| are small (scale-aware vs tol). */
ML_API ml_status_t ml_nelder_mead(ml_vec_func_t f, const double *x0, int n,
                                  double step, double tol, int max_iter,
                                  double *x_out, double *f_out);
/* ml_lbfgs_min: converges on ||g||<=tol*(1+||g0||); n<=32. */
ML_API ml_status_t ml_lbfgs_min(ml_vec_func_t f, const double *x0, int n,
                                double tol, int max_iter,
                                double *x_out, double *f_out);
/* ml_adam_min: full-batch finite-difference gradient per step, converges on
 * max|g|<=tol; bad lr returns ML_ERR_INVALID_ARG; n<=32. */
ML_API ml_status_t ml_adam_min(ml_vec_func_t f, const double *x0, int n,
                               double lr, double tol, int max_iter,
                               double *x_out, double *f_out);
#endif
