#ifndef LIBMATHC_QUADRATICS_H
#define LIBMATHC_QUADRATICS_H

#include "ml_core.h"

static inline double ml_equation(double a, double b, double c, double x) {
    return (a * x * x) + b * x + c;
}

// Citardauq Formula (Stable Quadratic) mapped to match naive +/- branches
static inline double ml_formula_pos(double a, double b, double c) {
    double disc = b * b - 4 * a * c;
    if (disc < 0.0) return 0.0 / 0.0;
    if (b >= 0.0) {
        double q = -0.5 * (b + ml_sqrt(disc));
        return c / q; // Maps to the +sqrt branch when b >= 0
    } else {
        double q = -0.5 * (b - ml_sqrt(disc));
        return q / a; // Maps to the +sqrt branch when b < 0
    }
}

static inline double ml_formula_neg(double a, double b, double c) {
    double disc = b * b - 4 * a * c;
    if (disc < 0.0) return 0.0 / 0.0;
    if (b >= 0.0) {
        double q = -0.5 * (b + ml_sqrt(disc));
        return q / a; // Maps to the -sqrt branch when b >= 0
    } else {
        double q = -0.5 * (b - ml_sqrt(disc));
        return c / q; // Maps to the -sqrt branch when b < 0
    }
}

// Legacy aliases for test suite compatibility
static inline double equation(double a, double b, double c, double x) { return ml_equation(a,b,c,x); }
static inline double formula_pos(double a, double b, double c) { return ml_formula_pos(a,b,c); }
static inline double formula_neg(double a, double b, double c) { return ml_formula_neg(a,b,c); }

#endif
