#ifndef MATHLIB_ML_EXP_LOG_H
#define MATHLIB_ML_EXP_LOG_H

#include "ml_core.h"

static inline double ml_exp(double x) {
    if (x == 0.0) return 1.0;
    if (x > 709.78) return 1.0 / 0.0;
    if (x < -745.13) return 0.0;
    double n = ml_round(x / 0.693147180559945309417);
    double r = x - n * 0.69314718036912381649 - n * 1.90821490974462528503e-10;
    // Horner's Method for exp(r) using precomputed inverse factorials
    // Eliminates division in the hot loop for maximum pipeline throughput
    static const double inv_fact[] = {
        1.0, 1.0, 0.5, 0.16666666666666666, 0.041666666666666664,
        0.008333333333333333, 0.001388888888888889, 0.0001984126984126984,
        2.48015873015873e-05, 2.7557319223985893e-06, 2.7557319223985888e-07,
        2.505210838544172e-08, 2.08767569878681e-09, 1.6059043836821613e-10,
        1.1470745597729725e-11, 7.647163731819816e-13, 4.779477332387385e-14,
        2.8114572543455206e-15, 1.5619206968586226e-16, 8.22063524662433e-18
    };
    double result = inv_fact[19];
    for (int i = 18; i >= 1; i--) result = __builtin_fma(result, r, inv_fact[i]);
    result = __builtin_fma(result, r, 1.0);
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

    // Horner's Method for 2 * atanh(z) extended to z^21 for < 1e-15 precision
    // Eliminates division and loop overhead for maximum pipeline throughput
    double poly = 0.09523809523809523; // 2/21
    poly = poly * z2 + 0.10526315789473684; // 2/19
    poly = poly * z2 + 0.11764705882352941; // 2/17
    poly = poly * z2 + 0.13333333333333333; // 2/15
    poly = poly * z2 + 0.15384615384615385; // 2/13
    poly = poly * z2 + 0.18181818181818182; // 2/11
    poly = poly * z2 + 0.2222222222222222;  // 2/9
    poly = poly * z2 + 0.2857142857142857;  // 2/7
    poly = poly * z2 + 0.4;                // 2/5
    poly = poly * z2 + 0.6666666666666666; // 2/3
    poly = poly * z2 + 2.0;                // 2/1

    return z * poly + e * 0.693147180559945309417;
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
