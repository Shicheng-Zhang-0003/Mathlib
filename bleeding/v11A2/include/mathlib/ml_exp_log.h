#ifndef MATHLIB_ML_EXP_LOG_H
#define MATHLIB_ML_EXP_LOG_H

#include "ml_core.h"

static inline double ml_exp(double x) {
    if (x == 0.0) return 1.0;
    if (x > 709.78) return 1.0 / 0.0;
    if (x < -745.13) return 0.0;

    double n = ml_round(x / 0.693147180559945309417);
    double r = x - n * 0.69314718036912381649 - n * 1.90821490974462528503e-10;

    double term = 1.0, result = 1.0;
    for (int i = 1; i <= 20; i++) {
        term *= r / i;
        result += term;
    }
    return ml_ldexp_pure(result, (int)n);
}

static inline double ml_log(double x) {
    if (x == 0.0) return -1.0 / 0.0;
    if (x < 0.0) return 0.0 / 0.0;
    if (x == 1.0) return 0.0;

    int e;
    double m = ml_frexp_pure(x, &e);
    if (m < 0.7071067811865475) { m *= 2.0; e--; }

    double z = (m - 1.0) / (m + 1.0);
    double z2 = z * z;
    double result = z, term = z;
    for (int i = 3; i <= 21; i += 2) {
        term *= z2;
        result += term / i;
    }
    return 2.0 * result + e * 0.693147180559945309417;
}

static inline double ml_pow(double x, double y) { return ml_exp(y * ml_log(x)); }


static inline double ml_logb(double x, double b) { return ml_log(x) / ml_log(b); }
static inline double ml_sinh(double x) { return (ml_exp(x) - ml_exp(-x)) / 2.0; }
static inline double ml_cosh(double x) { return (ml_exp(x) + ml_exp(-x)) / 2.0; }
static inline double ml_tanh(double x) { return ml_sinh(x) / ml_cosh(x); }
static inline double ml_asinh(double x) { return ml_log(x + ml_sqrt(x * x + 1.0)); }
static inline double ml_acosh(double x) { return (x < 1.0) ? 0.0/0.0 : ml_log(x + ml_sqrt(x * x - 1.0)); }
static inline double ml_atanh(double x) { return (x <= -1.0 || x >= 1.0) ? 0.0/0.0 : 0.5 * ml_log((1.0 + x) / (1.0 - x)); }

#endif
