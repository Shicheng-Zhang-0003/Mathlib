#ifndef LIBMATHC_INTEGRAL_H
#define LIBMATHC_INTEGRAL_H

#include "exponential.h"
#include "profiles.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_E
#define M_E 2.71828182845904523536
#endif
#define math_pi M_PI
#define math_e M_E

static inline double factorial_float(double x) {
    return (ml_sqrt(2 * math_pi * x)) * (power((x / math_e), x));
}

static inline double integral_traditional(double a, double b, double exponent, double additive, double d) {
    double result = 0.0;
    double x = a;
    while (x < b) {
        result += (power(x, exponent) + additive) * d;
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
