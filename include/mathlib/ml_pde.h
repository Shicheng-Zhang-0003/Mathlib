#ifndef MATHLIB_ML_PDE_H
#define MATHLIB_ML_PDE_H
#include "ml_compiler.h"
#include "ml_core.h"
#include "ml_types.h"
#ifdef __cplusplus
extern "C" {
#endif
/* One-step 1D PDE kernels: each call advances a SINGLE step (explicit heat,
 * implicit heat, wave) or solves one stationary problem (Poisson); the caller
 * loops for time integration. All boundary handling is Dirichlet-only:
 * endpoints are copied/held (u1[0]=u0[0], u1[n-1]=u0[n-1]; Poisson pins
 * u[0]=u[n-1]=0). No Neumann/Robin support. */
/* Explicit heat: one FTCS step. Requires r=alpha*dt/dx^2 <= 1/2 else
 * ML_ERR_SINGULAR. 3 <= n <= 4096. */
ML_API ml_status_t ml_heat_explicit(const double *u0, double *u1, int n, double dx, double dt, double alpha);
/* Implicit heat (backward Euler): one step via Thomas solve. Unstable-proof
 * for any r but still Dirichlet-only. 3 <= n <= 1024. */
ML_API ml_status_t ml_heat_implicit(const double *u0, double *u1, int n, double dx, double dt, double alpha);
/* Stationary Poisson -u''=f with Dirichlet 0 ends. 3 <= n <= 1024. */
ML_API ml_status_t ml_poisson_1d(const double *f, double *u, int n, double dx);
/* Second-order wave: one leapfrog step. Requires 0 <= cfl <= 1 else
 * ML_ERR_SINGULAR; Dirichlet ends held. 3 <= n <= 4096. */
ML_API ml_status_t ml_wave_leapfrog(const double *u_prev, const double *u_cur, double *u_next, int n, double cfl);
/* Single P1 element matrices (stiffness K_e, mass M_e) for mesh size h.
 * Element-local only: global assembly, Dirichlet BC application, and the
 * linear solve are the caller's responsibility. */
ML_API ml_status_t ml_fem1d_assemble(double h, double *ke00, double *ke01, double *ke11, double *me00, double *me01, double *me11);
#ifdef __cplusplus
}
#endif
#endif
