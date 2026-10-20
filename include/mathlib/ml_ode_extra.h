#ifndef MATHLIB_ML_ODE_EXTRA_H
#define MATHLIB_ML_ODE_EXTRA_H

#include "ml_compiler.h"
#include "ml_core.h"
#include "ml_types.h"
#include "ml_ode_sys.h"
#include "ml_numerical.h"
#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * v12R2 ODE / BVP DEPTH PACK
 *
 * Additional explicit and implicit integrators for first-order systems
 * (classic RK4, backward Euler, trapezoidal, BDF2, Adams-Bashforth and
 * Adams-Moulton predictor-corrector), a linear finite-difference boundary
 * value solver, and classic chaotic/oscillator right-hand sides.
 * ========================================================================== */

#define ML_ODE_MAXN 32

/* Explicit classic RK4 for y' = f(t,y).  Final state written to y_out. */
ML_API ml_status_t ml_ode_rk4_sys (ml_sys_func_t f, double t0,
                                   const double *y0, int n, double dt, int steps,
                                   void *ctx, double *y_out);

/* Implicit Euler with Newton iteration (numerical Jacobian). */
ML_API ml_status_t ml_ode_backward_euler_sys (ml_sys_func_t f, double t0,
                                              const double *y0, int n, double dt,
                                              int steps, void *ctx,
                                              double *y_out);

/* Implicit trapezoidal rule. */
ML_API ml_status_t ml_ode_trapezoidal_sys (ml_sys_func_t f, double t0,
                                           const double *y0, int n, double dt,
                                           int steps, void *ctx, double *y_out);

/* Implicit 2-step BDF (started with one backward-Euler step). */
ML_API ml_status_t ml_ode_bdf2_sys (ml_sys_func_t f, double t0,
                                    const double *y0, int n, double dt,
                                    int steps, void *ctx, double *y_out);

/* Explicit Adams-Bashforth of order 1..4 (started with RK4). */
ML_API ml_status_t ml_ode_adams_bashforth_sys (ml_sys_func_t f, double t0,
                                               const double *y0, int n,
                                               int order, double dt, int steps,
                                               void *ctx, double *y_out);

/* Adams-Bashforth-Moulton predictor-corrector (2nd order). */
ML_API ml_status_t ml_ode_adams_moulton_sys (ml_sys_func_t f, double t0,
                                             const double *y0, int n, double dt,
                                             int steps, void *ctx,
                                             double *y_out);

/* Linear BVP y'' + p(x) y' + q(x) y = r(x) on [x0,x1], Dirichlet y(x0)=ya,
 * y(x1)=yb.  x and y must hold n+2 entries (boundaries included). */
ML_API ml_status_t ml_bvp_fd_linear (ml_func_t p, ml_func_t q, ml_func_t r,
                                     double x0, double x1, double ya, double yb,
                                     int n, double *x, double *y);

/* Classic right-hand sides.  ctx may be NULL (defaults) or a double[3]. */
ML_API void ml_lorenz_rhs (double t, const double *y, double *dy, int n,
                           void *ctx);
ML_API void ml_rossler_rhs (double t, const double *y, double *dy, int n,
                            void *ctx);

#ifdef __cplusplus
}
#endif
#endif /* MATHLIB_ML_ODE_EXTRA_H */
