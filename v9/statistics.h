#ifndef LIBMATHC_STATISTICS_H
#define LIBMATHC_STATISTICS_H

//Library header file for statistics and probability
#include <math.h>
#include "combinatorics.h"
#include "exponential.h"
inline double mean (double *data, int n);
inline double mean (double *data, int n) {
    double sum = 0;
    for (int i = 0; i < n; i++) {sum += data[i];}
    return sum / n;
} inline double variance (double *data, int n);
inline double variance (double *data, int n) {
    double m = mean (data, n);
    double sum = 0;
    for (int i = 0; i < n; i++) {double diff = data[i] - m;
    sum += diff * diff;}
    return sum / n;
} inline double stddev (double *data, int n);
inline double stddev (double *data, int n) {return sqrt (variance (data, n));}
inline double binomial_pmf (int n, int k, double p);
inline double binomial_pmf (int n, int k, double p) {
    if (k < 0 || k > n || p < 0 || p > 1) {return 0.0 / 0.0;}
    return ncr (n, k) * pow (p, k) * pow (1 - p, n - k);
} inline double normal_pdf (double x, double mu, double sigma);
inline double normal_pdf (double x, double mu, double sigma) {
    if (sigma <= 0) {return 0.0 / 0.0;}
    double z = (x - mu) / sigma;
    return (1 / sqrt (2 * math_pi * sigma * sigma)) * exponential (-z * z / 2);
} inline void linear_regression (double *x, double *y, int n, double *out_m, double *out_b);
inline void linear_regression (double *x, double *y, int n, double *out_m, double *out_b) {
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
