#ifndef MATHLIB_ML_MCMC_H
#define MATHLIB_ML_MCMC_H
#include "ml_compiler.h"
#include "ml_core.h"
#include "ml_types.h"
typedef double (*ml_logpdf_t)(const double *x, int n, void *ctx);
/* Metropolis-Hastings random-walk sampler, n <= 16, second-half mean.
 * ctx is passed through to every logp call (logp(x, n, ctx)). */
ML_API ml_status_t ml_mh_sample_ctx(ml_logpdf_t logp, const double *x0, int n, uint64_t seed,
                                    double proposal_std, int steps, double *mean_out, double *acc_rate,
                                    void *ctx);
/* Legacy ABI wrapper: identical but passes NULL ctx to logp.
 * Prefer ml_mh_sample_ctx when the target needs user data. */
ML_API ml_status_t ml_mh_sample(ml_logpdf_t logp, const double *x0, int n, uint64_t seed,
                                double proposal_std, int steps, double *mean_out, double *acc_rate);
ML_API double ml_kde_gaussian(const double *data, int n, double x, double h);
ML_API double ml_ess(const double *xs, int n);
ML_API double ml_gelman_rubin(const double *chain1, const double *chain2, int n);
#endif
