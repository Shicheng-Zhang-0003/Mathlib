#ifndef MATHLIB_ML_STATS_INV_H
#define MATHLIB_ML_STATS_INV_H
#include "ml_compiler.h"
#include "ml_core.h"
ML_API double ml_gamma_inv(double p, double a, double b);
ML_API double ml_beta_inv(double p, double a, double b);
ML_API double ml_chi2_inv(double p, int k);
/* Starters: symmetric bracket around median 0, doubling expansion of the
 * open side until CDF brackets p (mirrors gamma/F inversion); <=300 CDF
 * evals total, then bisection. */
ML_API double ml_student_t_inv(double p, int nu);
ML_API double ml_f_inv(double p, int d1, int d2);
ML_API double ml_normal_logcdf(double x, double mu, double sigma);
#endif
