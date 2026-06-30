#ifndef LIBMATHC_EXPONENTIAL_H
#define LIBMATHC_EXPONENTIAL_H

#include "ml_core.h"
#ifndef M_E
#define M_E 2.71828182845904523536
#endif
#define math_e M_E
#define math_ln2 0.693147180559945309417

// Base-2 Split: e^x = 2^n * e^r
static inline double exponential(double x) {
    if (ml_isinf(x)) return (x > 0) ? x : 0.0;
    if (x == 0.0) return 1.0;
    double n = ml_round(x / math_ln2);
    double r = x - n * math_ln2;

    double result = 1.0;
    double term = 1.0;
    for (int i = 1; i <= 20; i++) {
        term *= r / i;
        result += term;
    }
    // ldexp efficiently multiplies by 2^n using IEEE 754 exponent bit-shifting
    return ldexp(result, (int)n);
}

// Fast Series: ln(x) = e * ln(2) + 2 * (z + z^3/3 + z^5/5...) where z = (m-1)/(m+1)
static inline double logarithm(double x) {
    if (x == 0.0) return -1.0 / 0.0;
    if (x < 0.0) return 0.0 / 0.0;
    if (x <= 0.0) return 0.0 / 0.0;
    if (x == 1.0) return 0.0;

    int e;
    // frexp extracts IEEE 754 exponent and mantissa (0.5 <= m < 1.0)
    double m = frexp(x, &e);

    // Adjust to [ml_sqrt(2)/2, ml_sqrt(2)] for optimal convergence
    if (m < 0.7071067811865475) {
        m *= 2.0;
        e -= 1;
    }

    double z = (m - 1.0) / (m + 1.0);
    double z2 = z * z;
    double result = z;
    double term = z;

    for (int i = 3; i <= 15; i += 2) {
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
