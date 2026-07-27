#ifndef MATHLIB_ML_STATISTICS_H
#define MATHLIB_ML_STATISTICS_H

#include "ml_compiler.h"
#include "ml_core.h"
#include "ml_combinatorics.h"
#include "ml_exp_log.h"

/* FIX: Added const to non-mutating data arrays */
ML_API double ml_mean(const double *data, int n);
ML_API double ml_variance(const double *data, int n);
ML_API double ml_stddev(const double *data, int n);
ML_API double ml_binomial_pmf(int n, int k, double p);
ML_API double ml_normal_pdf(double x, double mu, double sigma);
ML_API void ml_linear_regression(const double *x, const double *y, int n, double *out_m, double *out_b);

#endif /* MATHLIB_ML_STATISTICS_H */
