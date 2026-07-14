#ifndef LIBMATHC_INTEGRAL_H
#define LIBMATHC_INTEGRAL_H

#include "ml_exp_log.h"
#include "profiles.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_E
#define M_E 2.71828182845904523536
#endif
#ifndef math_pi
#define math_pi M_PI
#endif
#ifndef math_e
#define math_e M_E
#endif

static inline double factorial_float(double x) {
    return (ml_sqrt(2 * math_pi * x)) * (ml_pow((x / math_e), x));
}

static inline double integral_traditional(double a, double b, double exponent, double additive, double d) {
    // Temporal Termination: Hard ceiling to prevent infinite loops on malicious/NaN step sizes
    if (d == 0.0 || ml_isnan(d) || ml_isnan(a) || ml_isnan(b)) return 0.0/0.0;
    if ((d > 0 && a >= b) || (d < 0 && a <= b)) return 0.0;

    double result = 0.0;
    double x = a;
    int max_steps = 10000000; // O(1) temporal ceiling
    for (int step = 0; step < max_steps; step++) {
        if ((d > 0 && x >= b) || (d < 0 && x <= b)) break;
        result += (ml_pow(x, exponent) + additive) * d;
        x += d;
    }
    return result;
}

static inline double gamma_new(double x) {
    if (x <= 0) return 0.0 / 0.0;
    if (x > 2) return (x - 1) * gamma_new(x - 1);
    if (x < 1) return gamma_new(x + 1) / x;
    double z = x - 1;
    double p = -0.193527818 + z * 0.035868343;
    p = 0.482199394 + z * p;
    p = -0.756704078 + z * p;
    p = 0.918206857 + z * p;
    p = -0.897056937 + z * p;
    p = 0.989028236 + z * p;
    p = -0.577191652 + z * p;
    return 1 + z * p;
}
#endif
