#ifndef MATHLIB_ML_OPTIM_EXTRA_H
#define MATHLIB_ML_OPTIM_EXTRA_H

#include "ml_compiler.h"
#include "ml_core.h"
#include "ml_types.h"
#include "ml_optim_n.h"
#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * v12R2 OPTIMISATION DEPTH PACK
 *
 * Linear programming (primal simplex + branch & bound integer LP),
 * quadratic programming (equality KKT and primal active-set), constrained
 * nonlinear optimisation (augmented Lagrangian), Newton and Gauss-Newton,
 * and a projected-gradient box solver.
 * ========================================================================== */

/* ------------------------------------------------ linear programming ----- */
/* Maximise c^T x subject to A x <= b, x >= 0.  A is m x n row-major and b
 * must be nonnegative.  x receives an optimal vertex, *obj the optimum.
 * Unbounded problems return ML_ERR_CONVERGENCE. */
ML_API ml_status_t ml_simplex_lp (const double *c, const double *A,
                                  const double *b, int m, int n, double *x,
                                  double *obj);

/* Minimise c^T x subject to A x <= b, x >= 0. */
ML_API ml_status_t ml_simplex_lp_min (const double *c, const double *A,
                                      const double *b, int m, int n, double *x,
                                      double *obj);

/* Integer LP: maximise c^T x, A x <= b, x >= 0, x integer. */
ML_API ml_status_t ml_integer_lp (const double *c, const double *A,
                                  const double *b, int m, int n, double *x,
                                  double *obj, int max_nodes);

/* ------------------------------------------------- quadratic programming - */
/* Minimise 1/2 x^T H x + g^T x subject to A x = b (A is m x n). */
ML_API ml_status_t ml_qp_eq (const double *H, const double *g, const double *A,
                             const double *b, int n, int m, double *x);

/* Minimise 1/2 x^T H x + g^T x subject to A x <= b (A is m x n). */
ML_API ml_status_t ml_qp_active_set (const double *H, const double *g,
                                     const double *A, const double *b, int n,
                                     int m, const double *x0, int max_iter,
                                     double tol, double *x);

/* ----------------------------------------- constrained nonlinear --------- */
/* Residual callback for nonlinear least squares: fill r[0..m-1]. */
typedef void (*ml_residual_t) (const double *x, int n, double *r, void *ctx);

/* Gauss-Newton / Levenberg-Marquardt on sum r_i(x)^2. */
ML_API ml_status_t ml_gauss_newton (ml_residual_t res, int m, const double *x0,
                                    int n, int max_iter, double tol,
                                    double *x_out);

/* Unconstrained Newton with finite-difference gradient/Hessian. */
ML_API ml_status_t ml_newton_min (ml_vec_func_t f, const double *x0, int n,
                                  int max_iter, double tol, double *x_out,
                                  double *f_out);

/* Box-constrained minimisation low<=x<=up (ML_GRAPH-free projected gradient).
 * Pass +-inf to leave a bound inactive. */
ML_API ml_status_t ml_projected_gradient_box (ml_vec_func_t f,
                                              const double *x0, const double *lo,
                                              const double *up, int n,
                                              int max_iter, double tol,
                                              double *x_out, double *f_out);

/* Augmented Lagrangian for  min f(x)  s.t.  h_i(x)=0,  gj(x)<=0.
 * h and g may be NULL when nh/ng are zero. */
ML_API ml_status_t ml_augmented_lagrangian (ml_vec_func_t f, ml_residual_t h,
                                            int nh, ml_residual_t g, int ng,
                                            void *hctx, void *gctx,
                                            const double *x0, int n,
                                            int max_iter, double tol,
                                            double *x_out, double *f_out);

#ifdef __cplusplus
}
#endif
#endif /* MATHLIB_ML_OPTIM_EXTRA_H */
