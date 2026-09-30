#ifndef MATHLIB_ML_STATISTICS_H
#define MATHLIB_ML_STATISTICS_H

#include "ml_compiler.h"
#include "ml_core.h"
#include "ml_combinatorics.h"
#include "ml_exp_log.h"
#ifdef __cplusplus
extern "C" {
#endif

/* FIX: Added const to non-mutating data arrays */
ML_API double ml_mean(const double *data, int n);
/* Population variance (divide by n). Sample variant below divides by n-1;
 * ml_covariance is sample (n-1). ml_stddev wraps population variance. */
ML_API double ml_variance(const double *data, int n);
ML_API double ml_variance_s(const double *data, int n);
ML_API double ml_stddev(const double *data, int n);
/* tmp must hold at least n doubles (scratch sort buffer, tmp[n]). */
ML_API double ml_median(const double *data, int n, double *tmp);
ML_API double ml_covariance(const double *x, const double *y, int n);
ML_API double ml_correlation(const double *x, const double *y, int n);
ML_API double ml_hmean(const double *data, int n);
ML_API double ml_gmean(const double *data, int n);
ML_API double ml_geometric_pmf(int k, double p);
ML_API double ml_uniform_pdf(double x, double a, double b);
ML_API double ml_binom_mean(int n, double p);
ML_API double ml_binom_var(int n, double p);
ML_API double ml_bayes(double p_b_given_a, double p_a, double p_b);
ML_API double ml_binomial_pmf(int n, int k, double p);
ML_API double ml_normal_pdf(double x, double mu, double sigma);
ML_API double ml_normal_cdf(double x, double mu, double sigma);
ML_API double ml_normal_inv(double p, double mu, double sigma);
ML_API double ml_exponential_pdf(double x, double lambda);
ML_API double ml_poisson_pmf(int k, double lambda);
ML_API void ml_linear_regression(const double *x, const double *y, int n, double *out_m, double *out_b);
/* Undergraduate probability: exact transforms of gamma_p/q and beta.
 * chi2(k,x)=P(k/2,x/2); exponential_cdf=1-e^{-lx}; gamma_cdf=P(a,bx);
 * beta/Student-t/F via hyp2f1-free continued fractions in integral.c. */
ML_API double ml_gamma_cdf(double x, double a, double b);
ML_API double ml_chi2_cdf(double x, int k);
ML_API double ml_exponential_cdf(double x, double lambda);
ML_API double ml_beta_cdf(double x, double a, double b);
ML_API double ml_student_t_cdf(double t, int nu);
ML_API double ml_f_cdf(double x, int d1, int d2);

#ifdef __cplusplus
}
#endif
#endif /* MATHLIB_ML_STATISTICS_H */
