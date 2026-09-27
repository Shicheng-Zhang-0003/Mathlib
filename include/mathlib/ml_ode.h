#ifndef MATHLIB_ML_ODE_H
#define MATHLIB_ML_ODE_H

#include "ml_compiler.h"
#include "ml_core.h"

typedef double (*ml_ode_func_t)(double t, double y);

ML_API double ml_ode_euler(ml_ode_func_t f, double t0, double y0, double dt, int steps);
ML_API double ml_ode_rk4(ml_ode_func_t f, double t0, double y0, double dt, int steps);
ML_API double ml_ode_dp5(ml_ode_func_t f, double t0, double y0, double t1, double rtol, double atol);
/* Explicit midpoint RK2: y_{n+1} = y_n + dt*f(t+dt/2, y+dt/2*f(t,y)).
 * Second-order but NOT symplectic and NOT energy-preserving (it is not
 * velocity-Verlet / leapfrog in the Hamiltonian sense, which needs a
 * position/velocity split). Do not use for energy conservation. */
ML_API double ml_ode_midpoint2(ml_ode_func_t f, double t0, double y0, double dt, int steps);
/* DEPRECATED alias of ml_ode_midpoint2, kept for ABI. The old name is a
 * misnomer: this is explicit midpoint, not symplectic, do not use for energy. */
ML_API ML_DEPRECATED("ml_ode_leapfrog is explicit midpoint, not symplectic, do not use for energy; use ml_ode_midpoint2")
double ml_ode_leapfrog(ml_ode_func_t f, double t0, double y0, double dt, int steps);
ML_API double ml_ode_heun(ml_ode_func_t f, double t0, double y0, double dt, int steps);
ML_API double ml_suvat_s(double u, double a, double t);
ML_API double ml_suvat_v(double u, double a, double t);
ML_API double ml_projectile_range(double v0, double theta, double g);
ML_API double ml_projectile_height(double v0, double theta, double g);
ML_API double ml_shm_x(double A, double omega, double t, double phi);
ML_API double ml_collision_1d(double m1, double m2, double u1, double u2, double *v1, double *v2);
ML_API int ml_ode2_const(double a, double b, double c, double *r0, double *r1);

#endif /* MATHLIB_ML_ODE_H */
