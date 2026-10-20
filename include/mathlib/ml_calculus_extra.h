#ifndef MATHLIB_ML_CALCULUS_EXTRA_H
#define MATHLIB_ML_CALCULUS_EXTRA_H

#include "ml_compiler.h"
#include "ml_core.h"
#include "ml_types.h"
#include "ml_numerical.h"
#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * v12R2 CALCULUS DEPTH PACK
 *
 * Multivariate finite differences, limits, series acceleration, Pade
 * approximants and multidimensional / line / surface quadrature on top of the
 * scalar routines in ml_numerical.h and ml_calculus.h.
 * ========================================================================== */

/* R^n -> R and R^n -> R^m callbacks with an opaque context. */
typedef double (*ml_rn_func_t) (const double *x, int n, void *ctx);
typedef void (*ml_rn_vec_func_t) (const double *x, int n, double *out, int m,
                                  void *ctx);
typedef double (*ml_r2_func_t) (double x, double y, void *ctx);
typedef double (*ml_r1c_func_t) (double x, void *ctx);
typedef double (*ml_r3c_func_t) (double x, double y, double z, void *ctx);
typedef double (*ml_uv_func_t) (double u, double v);

/* Richardson-extrapolated first derivative of a scalar function.
 * *err_out (optional) receives an error estimate. */
ML_API double ml_derivative_richardson (ml_func_t f, double x, double h,
                                        double *err_out);

/* Gradient (length n), Jacobian (m x n, row-major) and Hessian (n x n) by
 * central differences / Richardson. All require finite x. */
ML_API ml_status_t ml_gradient_rn (ml_rn_func_t f, const double *x, int n,
                                   double *grad, double h, void *ctx);
ML_API ml_status_t ml_jacobian_rn (ml_rn_vec_func_t f, const double *x, int n,
                                   int m, double *J, double h, void *ctx);
ML_API ml_status_t ml_hessian_rn (ml_rn_func_t f, const double *x, int n,
                                  double *H, double h, void *ctx);
ML_API double ml_directional_derivative (ml_rn_func_t f, const double *x, int n,
                                         const double *dir, double h, void *ctx);

/* One-sided (side = +1 right, -1 left) or two-sided (side = 0) limit of
 * f as x -> x0, evaluated by Richardson collapse. *err_out optional. */
ML_API double ml_limit (ml_func_t f, double x0, int side, double tol,
                        double *err_out);

/* Series acceleration. ml_aitken applies one Aitken Delta^2 step to the last
 * three partial sums. ml_richardson_geom assumes partial sums obey
 * s_n = s + C q^n and eliminates the dominant error. */
ML_API double ml_aitken (const double *s, int n);
ML_API double ml_richardson_geom (const double *s, int n, double q);

/* Pade approximant [L/M] of a series with Taylor coefficients c[0..L+M].
 * p has L+1 entries, q has M+1 entries with q[0] = 1. */
ML_API ml_status_t ml_pade (const double *c, int L, int M, double *p, double *q);
ML_API double ml_pade_eval (const double *p, const double *q, int L, int M,
                            double x);

/* Tensor-product composite Simpson quadrature. n counts intervals (even). */
ML_API double ml_integral_2d (ml_r2_func_t f, double xa, double xb, double ya,
                              double yb, int nx, int ny, void *ctx);
ML_API double ml_integral_2d_var (ml_r2_func_t f, double xa, double xb,
                                  ml_r1c_func_t ylo, ml_r1c_func_t yhi,
                                  int nx, int ny, void *ctx);
ML_API double ml_integral_3d (ml_r3c_func_t f, double xa, double xb, double ya,
                              double yb, double za, double zb, int n, void *ctx);

/* Line integrals of a scalar field and of a vector field along the curve
 * r(t) = (rx(t), ry(t), rz(t)), t in [a,b]. */
ML_API double ml_arc_length (ml_func_t rx, ml_func_t ry, ml_func_t rz,
                             double a, double b, int n);
ML_API double ml_line_integral_scalar (ml_r3c_func_t f, ml_func_t rx,
                                       ml_func_t ry, ml_func_t rz, double a,
                                       double b, int n, void *ctx);
ML_API double ml_line_integral_vector (ml_r3c_func_t fx, ml_r3c_func_t fy,
                                       ml_r3c_func_t fz, ml_func_t rx,
                                       ml_func_t ry, ml_func_t rz, double a,
                                       double b, int n, void *ctx);

/* Surface integral of a scalar field over r(u,v), (u,v) in the rectangle. */
ML_API double ml_surface_integral_scalar (ml_r3c_func_t f, ml_uv_func_t rx,
                                          ml_uv_func_t ry, ml_uv_func_t rz,
                                          double ua, double ub, double va,
                                          double vb, int nu, int nv, void *ctx);

#ifdef __cplusplus
}
#endif
#endif /* MATHLIB_ML_CALCULUS_EXTRA_H */
