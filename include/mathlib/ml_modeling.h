#ifndef MATHLIB_ML_MODELING_H
#define MATHLIB_ML_MODELING_H

#include "ml_compiler.h"
#include "ml_core.h"
#include "ml_types.h"
#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Applied modeling kit: Markov chains, epidemics, populations, queues,
 * conjugate Bayes, chaos. The postgrad-contest (COMAP-style) companion:
 * closed forms where they exist, small verified solvers elsewhere.
 *
 * Conventions: row-stochastic matrices (rows sum to 1); rates >= 0 with
 * stability required (lambda < mu for queues, R0 documented for SIR);
 * NaN on invalid input; iterative solvers return ML_ERR_CONVERGENCE on
 * exhaustion (never silent SUCCESS). Caps keep every frame < 8KB.
 * ========================================================================== */

/* Markov chains (P row-stochastic, n <= 256) */
ML_API int ml_markov_is_stochastic(const double *P, int n, double tol);
ML_API ml_status_t ml_markov_step(const double *P, const double *pi,
                                  double *out, int n);
ML_API ml_status_t ml_markov_stationary(const double *P, double *pi, int n,
                                        int max_iter, double tol);

/* SIR epidemiology (fractions of population) */
ML_API double ml_sir_r0(double beta, double gamma);
ML_API double ml_sir_final_size(double r0);
ML_API double ml_sir_peak(double r0);

/* Leslie population (first row fecundity >= 0, subdiagonal survival [0,1]) */
ML_API ml_status_t ml_leslie_step(const double *L, const double *pop,
                                  double *out, int n);
ML_API ml_status_t ml_leslie_growth(const double *L, int n, double *lambda);

/* M/M/1 queues (0 <= lambda < mu) */
ML_API double ml_mm1_rho(double lambda, double mu);
ML_API double ml_mm1_L(double lambda, double mu);
ML_API double ml_mm1_W(double lambda, double mu);
ML_API double ml_mm1_wait(double lambda, double mu, double t);

/* Conjugate Bayes (exact u64 counts; intervals via staged beta_inv) */
ML_API ml_status_t ml_beta_binomial_post(uint64_t a0, uint64_t b0,
                                         uint64_t succ, uint64_t trials,
                                         double *a1, double *b1);
ML_API double ml_beta_post_mean(double a1, double b1);
ML_API ml_status_t ml_gamma_poisson_post(double a0, double b0, uint64_t k,
                                         double texp, double *a1, double *b1);
ML_API ml_status_t ml_beta_credible(double a, double b, double alpha,
                                    double *lo, double *hi);

/* Chaos: logistic-map Lyapunov exponent estimate (r in [0,4]) */
ML_API double ml_logistic_lyapunov(double r, double x0, int transient, int steps);

#ifdef __cplusplus
}
#endif
#endif /* MATHLIB_ML_MODELING_H */
