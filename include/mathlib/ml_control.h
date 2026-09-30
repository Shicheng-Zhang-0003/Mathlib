#ifndef MATHLIB_ML_CONTROL_H
#define MATHLIB_ML_CONTROL_H
#include "ml_compiler.h"
#include "ml_core.h"
#include "ml_types.h"
#include "ml_tensor.h"
#ifdef __cplusplus
extern "C" {
#endif
/* Continuous-time LQR for the 2-state SISO plant dx=(A)x+(B)u with scalar
 * weights Q=q*I (q>0) and R=r (r>0): NO full Q matrix, NO cross-weight N,
 * SISO-only (single input). Solves the ARE via Newton-Kleinman; returns
 * ML_ERR_SINGULAR on singular Lyapunov steps AND on iteration exhaustion
 * (no dedicated non-convergence code; never returns success un-converged). */
ML_API ml_status_t ml_lqr_gain_2x2(double a00, double a01, double a10, double a11,
                                   double b0, double b1, double q, double r,
                                   double *k0, double *k1);
ML_API ml_status_t ml_kalman_1d(double x0, double p0, const double *zs, int n,
                                double F, double H, double Q, double R,
                                double *x_out, double *p_out);
ML_API double ml_lyapunov_2x2_trace(double a00, double a01, double a10, double a11);
#ifdef __cplusplus
}
#endif
#endif
