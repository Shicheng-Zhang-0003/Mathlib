#ifndef LIBMATHC_STATISTICS_H
#define LIBMATHC_STATISTICS_H

//Library header file for statistics and probability
#include "ml_core.h"
#include "combinatorics.h"
#include "ml_exp_log.h"
static inline double mean (double *data, int n);
static inline double mean (double *data, int n) {
    double sum = 0;
    for (int i = 0; i < n; i++) {sum += data[i];}
    return sum / n;
} static inline double variance (double *data, int n);
static inline double variance (double *data, int n) {
    double m = mean (data, n);
    double sum = 0;
    for (int i = 0; i < n; i++) {double diff = data[i] - m;
    sum += diff * diff;}
    return sum / n;
} static inline double stddev (double *data, int n);
static inline double stddev (double *data, int n) {return ml_sqrt(variance (data, n));}
static inline double binomial_pmf (int n, int k, double p);
static inline double binomial_pmf (int n, int k, double p) {
    if (k < 0 || k > n || p < 0 || p > 1) {return 0.0 / 0.0;}
    return ncr (n, k) * ml_pow(p, k) * ml_pow(1 - p, n - k);
} static inline double normal_pdf (double x, double mu, double sigma);
static inline double normal_pdf (double x, double mu, double sigma) {
    if (sigma <= 0) {return 0.0 / 0.0;}
    double z = (x - mu) / sigma;
    return (1 / ml_sqrt(2 * math_pi * sigma * sigma)) * ml_exp (-z * z / 2);
} static inline void linear_regression (double *x, double *y, int n, double *out_m, double *out_b);
static inline void linear_regression (double *x, double *y, int n, double *out_m, double *out_b) {
    double sum_x = 0, sum_y = 0, sum_xy = 0, sum_x2 = 0;
    for (int i = 0; i < n; i++) {sum_x += x[i];
    sum_y += y[i];
    sum_xy += x[i] * y[i];
    sum_x2 += x[i] * x[i];}
    double denom = n * sum_x2 - sum_x * sum_x;
    if (denom == 0) {*out_m = 0.0 / 0.0; *out_b = 0.0 / 0.0; return;}
    *out_m = (n * sum_xy - sum_x * sum_y) / denom;
    *out_b = (sum_y * sum_x2 - sum_x * sum_xy) / denom;
}



#endif
