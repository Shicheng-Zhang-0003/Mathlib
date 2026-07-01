#ifndef LIBMATHC_EXPONENTIAL_H
#define LIBMATHC_EXPONENTIAL_H

#include "ml_core.h"

#define math_ln2 0.693147180559945309417
#define LN2_HI 0.69314718036912381649017333984375
#define LN2_LO 1.9082149097446252850341796875e-10

static inline double exponential(double x) {
    if (x == 0.0) return 1.0;
    if (x > 709.78) return 1.0 / 0.0;
    if (x < -745.13) return 0.0;

    double n = ml_round(x / math_ln2);
    double r = x - n * LN2_HI - n * LN2_LO; // Cody-Waite split

    double term = 1.0;
    double result = 1.0;
    for (int i = 1; i <= 20; i++) {
        term *= r / i;
        result += term;
    }
    return ml_ldexp_pure(result, (int)n);
}

static inline double logarithm(double x) {
    if (x == 0.0) return -1.0 / 0.0;
    if (x < 0.0) return 0.0 / 0.0;
    if (x == 1.0) return 0.0;

    int e;
    double m = ml_frexp_pure(x, &e);

    // CRITICAL FIX: Adjust m to [sqrt(2)/2, sqrt(2)] to minimize z
    if (m < 0.7071067811865475) {
        m *= 2.0;
        e--;
    }

    double z = (m - 1.0) / (m + 1.0);
    double z2 = z * z;
    double result = z;
    double term = z;
    for (int i = 3; i <= 21; i += 2) {
        term *= z2;
        result += term / i;
    }
    return 2.0 * result + e * math_ln2;
}

static inline double power(double x, double y) { return exponential(y * logarithm(x)); }
static inline double logarithm_base(double x, double b) { return logarithm(x) / logarithm(b); }
static inline double hyperbolic_sine(double x) { return (exponential(x) - exponential(-x)) / 2.0; }
static inline double hyperbolic_cosine(double x) { return (exponential(x) + exponential(-x)) / 2.0; }
static inline double hyperbolic_tangent(double x) { return hyperbolic_sine(x) / hyperbolic_cosine(x); }
static inline double inverse_hyperbolic_sine(double x) { return logarithm(x + ml_sqrt(x * x + 1.0)); }
static inline double inverse_hyperbolic_cosine(double x) { return (x < 1.0) ? 0.0/0.0 : logarithm(x + ml_sqrt(x * x - 1.0)); }
static inline double inverse_hyperbolic_tangent(double x) { return (x <= -1.0 || x >= 1.0) ? 0.0/0.0 : 0.5 * logarithm((1.0 + x) / (1.0 - x)); }

#endif
