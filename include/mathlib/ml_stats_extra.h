#ifndef MATHLIB_ML_STATS_EXTRA_H
#define MATHLIB_ML_STATS_EXTRA_H

#include "ml_compiler.h"
#include "ml_core.h"
#include "ml_types.h"
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * v12R2 STATISTICS DEPTH PACK
 *
 * Descriptive statistics, additional distributions, hypothesis tests with
 * p-values, confidence intervals, (multiple / ridge / logistic) regression,
 * maximum-likelihood estimation and bootstrap intervals.
 * Builds on the cdf/inverse machinery in ml_statistics.h / ml_stats_inv.h.
 * ========================================================================== */

/* ------------------------------------------------------- descriptive ----- */
ML_API double ml_stat_skewness (const double *d, int n);
ML_API double ml_stat_kurtosis (const double *d, int n); /* excess */
ML_API double ml_quantile (const double *d, int n, double q, double *tmp);
ML_API double ml_percentile (const double *d, int n, double p, double *tmp);
ML_API double ml_iqr (const double *d, int n, double *tmp);
ML_API double ml_cv (const double *d, int n);
ML_API double ml_range (const double *d, int n);
ML_API double ml_mad (const double *d, int n, double *tmp);
ML_API double ml_trimmed_mean (const double *d, int n, double alpha, double *tmp);
ML_API ml_status_t ml_cov_matrix (const double *const *cols, int n, int p,
                                  double *cov);
ML_API ml_status_t ml_corr_matrix (const double *const *cols, int n, int p,
                                   double *corr);

/* ------------------------------------------------------ distributions ---- */
ML_API double ml_lognormal_pdf (double x, double mu, double sigma);
ML_API double ml_lognormal_cdf (double x, double mu, double sigma);
ML_API double ml_lognormal_inv (double p, double mu, double sigma);
ML_API double ml_weibull_pdf (double x, double shape, double scale);
ML_API double ml_weibull_cdf (double x, double shape, double scale);
ML_API double ml_cauchy_pdf (double x, double loc, double scale);
ML_API double ml_cauchy_cdf (double x, double loc, double scale);
ML_API double ml_uniform_cdf (double x, double a, double b);
ML_API double ml_beta_pdf (double x, double a, double b);
ML_API double ml_gamma_pdf (double x, double shape, double rate);
ML_API double ml_chi2_pdf (double x, double k);
ML_API double ml_student_t_pdf (double x, double nu);
ML_API double ml_f_pdf (double x, double d1, double d2);
ML_API double ml_poisson_cdf (int k, double lambda);
ML_API double ml_binomial_cdf (int n, int k, double p);
ML_API double ml_bernoulli_pmf (int k, double p);
ML_API double ml_negbinom_pmf (int k, double r, double p);
ML_API double ml_hypergeometric_pmf (int k, int N, int K, int n);
ML_API double ml_multivariate_normal_pdf (const double *x, const double *mu,
                                          const double *Sigma, int n,
                                          double *work);

/* -------------------------------------------------- hypothesis tests ----- */
ML_API ml_status_t ml_ttest_1samp (const double *d, int n, double mu0,
                                   double *t, double *p);
ML_API ml_status_t ml_ttest_2samp (const double *x, int nx, const double *y,
                                   int ny, double *t, double *p);
ML_API ml_status_t ml_ttest_welch (const double *x, int nx, const double *y,
                                   int ny, double *t, double *p);
ML_API ml_status_t ml_ttest_paired (const double *x, const double *y, int n,
                                    double *t, double *p);
ML_API ml_status_t ml_ztest_mean (const double *d, int n, double mu0,
                                  double sigma, double *z, double *p);
ML_API ml_status_t ml_f_test_var (const double *x, int nx, const double *y,
                                  int ny, double *F, double *p);
ML_API ml_status_t ml_chi2_gof (const int *obs, const double *exp, int k,
                                double *chi2, double *p);
ML_API ml_status_t ml_anova_oneway (const double *const *groups,
                                    const int *sizes, int g, double *F,
                                    double *p);
ML_API ml_status_t ml_mann_whitney_u (const double *x, int nx, const double *y,
                                      int ny, double *U, double *p);

/* ------------------------------------------------- confidence intervals - */
ML_API ml_status_t ml_ci_mean (const double *d, int n, double level, double *lo,
                               double *hi);
ML_API ml_status_t ml_ci_proportion (int k, int n, double level, double *lo,
                                     double *hi);
ML_API ml_status_t ml_ci_variance (const double *d, int n, double level,
                                   double *lo, double *hi);

/* ------------------------------------------------------- regression ------ */
/* X is n x p row-major (caller includes an intercept column if desired). */
ML_API ml_status_t ml_multiple_regression (const double *X, const double *y,
                                           int n, int p, double *beta,
                                           double *sigma2);
ML_API ml_status_t ml_ridge_regression (const double *X, const double *y, int n,
                                        int p, double lambda, double *beta);
ML_API ml_status_t ml_wls (const double *X, const double *y, const double *w,
                           int n, int p, double *beta);
ML_API ml_status_t ml_logistic_regression (const double *X, const double *y,
                                           int n, int p, double *beta,
                                           int max_iter, double tol);

/* ------------------------------------------------------------- MLE ------- */
ML_API ml_status_t ml_mle_normal (const double *d, int n, double *mu,
                                  double *sigma2);
ML_API ml_status_t ml_mle_exponential (const double *d, int n, double *lambda);
ML_API ml_status_t ml_mle_poisson (const int *d, int n, double *lambda);
ML_API ml_status_t ml_mle_bernoulli (const int *d, int n, double *p);

/* --------------------------------------------------------- bootstrap ----- */
ML_API ml_status_t ml_bootstrap_ci_mean (const double *d, int n, int B,
                                         double level, uint64_t seed,
                                         double *lo, double *hi);

#ifdef __cplusplus
}
#endif
#endif /* MATHLIB_ML_STATS_EXTRA_H */
