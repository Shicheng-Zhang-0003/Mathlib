#ifndef MATHLIB_ML_SDE_H
#define MATHLIB_ML_SDE_H
#include "ml_compiler.h"
#include "ml_core.h"
#include "ml_types.h"
typedef double (*ml_sde_drift_t)(double t, double x, void *ctx);
typedef double (*ml_sde_diff_t)(double t, double x, void *ctx);
ML_API ml_status_t ml_sde_euler_maruyama(ml_sde_drift_t a, ml_sde_diff_t b, double t0, double x0,
                                         double t1, int steps, uint64_t seed, double *x1);
ML_API ml_status_t ml_sde_milstein(ml_sde_drift_t a, ml_sde_diff_t b,
                                   double (*dbdx)(double t, double x, void *ctx),
                                   double t0, double x0, double t1, int steps, uint64_t seed, double *x1);
/* Mean-only pinned bridge E[W_t | W_0=a, W_T=wT] = a+(wT-a)*t/T.
 * b is kept for API symmetry and ignored (must equal wT for a pinned
 * bridge). Deterministic: for a random draw use ml_brownian_bridge_sample. */
ML_API double ml_brownian_bridge(double t, double T, double a, double b, double wT);
/* Honest alias of ml_brownian_bridge (mean-only, no variance term). */
ML_API double ml_brownian_bridge_mean(double t, double T, double a, double b, double wT);
/* TRUE Brownian bridge sampler: W_t | W_0=a, W_T=b ~
 * N(a+(b-a)*t/T, t*(T-t)/T), Gaussian from the module xorshift+Box-Muller
 * stream, deterministic in seed (seed 0 falls back to a fixed nonzero). */
ML_API double ml_brownian_bridge_sample(double a, double b, double T, double t, uint64_t seed);
ML_API double ml_ou_exact(double x0, double theta, double mu, double sigma, double t, double dw);
#endif
