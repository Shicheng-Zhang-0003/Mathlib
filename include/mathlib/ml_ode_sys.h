#ifndef MATHLIB_ML_ODE_SYS_H
#define MATHLIB_ML_ODE_SYS_H
#include "ml_compiler.h"
#include "ml_core.h"
#include "ml_types.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef void (*ml_sys_func_t)(double t, const double *y, double *dy, int n, void *ctx);
/* Explicit/DP5 system solver. n <= 16. */
ML_API ml_status_t ml_ode_dp5_sys(ml_sys_func_t f, double t0, const double *y0, int n,
                                  double t1, double rtol, double atol, double *y1, void *ctx);
/* Backward-Euler step-doubling (L-stable, not a ROS23 tableau). n <= 8.
 * One backward-Euler Newton solve per full/half step with finite-difference
 * Jacobian; error estimate from full-step vs two-half-steps difference. */
ML_API ml_status_t ml_ode_be2_sys(ml_sys_func_t f, double t0, const double *y0, int n,
                                  double t1, double rtol, double atol, double *y1, void *ctx);
/* DEPRECATED alias of ml_ode_be2_sys, kept for ABI. The old name is a
 * misnomer: backward-Euler step-doubling (L-stable, not ROS23 tableau). */
ML_API ML_DEPRECATED("ml_ode_ros23_sys is backward-Euler step-doubling, not a ROS23 tableau; use ml_ode_be2_sys")
ml_status_t ml_ode_ros23_sys(ml_sys_func_t f, double t0, const double *y0, int n,
                             double t1, double rtol, double atol, double *y1, void *ctx);
/* True symplectic (velocity-Verlet) splitting for separable Hamiltonians.
 * n <= 16. This, not ml_ode_midpoint2, is the energy-behaved integrator. */
ML_API ml_status_t ml_ode_symplectic_verlet(void (*acc)(const double *q, double *a, int n, void *ctx),
                                            double *q, double *p, int n, double h, int steps, void *ctx);
#ifdef __cplusplus
}
#endif
#endif
