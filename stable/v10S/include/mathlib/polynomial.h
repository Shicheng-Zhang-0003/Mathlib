#ifndef LIBMATHC_POLYNOMIAL_H
#define LIBMATHC_POLYNOMIAL_H

//Library header file for polynomial operations
#include "ml_core.h"
static inline double polynomial_eval (double *coeffs, int degree, double x);
static inline double polynomial_eval (double *coeffs, int degree, double x) {
    double result = coeffs[degree];
    for (int i = degree - 1; i >= 0; i--) {result = result * x + coeffs[i];}
    return result;
} static inline void polynomial_derivative (double *coeffs, int degree, double *out);
static inline void polynomial_derivative (double *coeffs, int degree, double *out) {
    for (int i = 0; i < degree; i++) {out[i] = coeffs[i + 1] * (i + 1);}
} static inline double polynomial_newton (double *coeffs, int degree, double x0, double epsilon, int max_iter);
static inline double polynomial_newton (double *coeffs, int degree, double x0, double epsilon, int max_iter) {
    double x = x0;
    for (int iter = 0; iter < max_iter; iter++) {
        double fx = coeffs[degree];
        double dfx = 0;
        for (int i = degree - 1; i >= 0; i--) {dfx = dfx * x + fx;
        fx = fx * x + coeffs[i];}
        if (ml_fabs(dfx) < epsilon) {return 0.0 / 0.0;}
        double x_next = x - fx / dfx;
        if (ml_fabs(x_next - x) < epsilon) {return x_next;}
        x = x_next;
    } return x;
}



#endif
