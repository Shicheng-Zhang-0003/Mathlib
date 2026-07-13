#ifndef MATHLIB_ML_TRIG_H
#define MATHLIB_ML_TRIG_H

#include "ml_core.h"
#include "internal/minimax.h"
#include "internal/cordic.h"

#ifndef math_pi
#define math_pi 3.14159265358979323846
#endif

// Profile-Routed Trigonometry (The v11A2 Way)
static inline double ml_sin(double x) {
    if (ml_isnan(x) || ml_isinf(x)) return 0.0/0.0;
#if defined(MATHLIB_PROFILE_EMBEDDED)
    double s, c;
    ml_cordic_sincos(x, &s, &c);
    return s;
#else
    return ml_minimax_sin(x);
#endif
}

static inline double ml_cos(double x) {
    if (ml_isnan(x) || ml_isinf(x)) return 0.0/0.0;
#if defined(MATHLIB_PROFILE_EMBEDDED)
    double s, c;
    ml_cordic_sincos(x, &s, &c);
    return c;
#else
    return ml_minimax_cos(x);
#endif
}

static inline double ml_tan(double x) {
    double s = ml_sin(x);
    double c = ml_cos(x);
    if (c == 0.0) return (s > 0) ? 1.0/0.0 : -1.0/0.0;
    return s / c;
}

// --- Inverse Trigonometry ---

static inline double ml_atan(double x) {
    if (x > 1.0) return (math_pi / 2.0) - ml_atan(1.0 / x);
    if (x < -1.0) return -(math_pi / 2.0) - ml_atan(1.0 / x);
    if (x > 0.5) return (math_pi / 4.0) + ml_atan((x - 1.0) / (x + 1.0));
    if (x < -0.5) return -(math_pi / 4.0) + ml_atan((x + 1.0) / (1.0 - x));
    double result = x, term = x, x2 = x * x;
    for (int i = 3; i <= 21; i += 2) { term *= -x2; result += term / i; }
    return result;
}

static inline double ml_asin(double x) {
    if (x < -1.0 || x > 1.0) return 0.0 / 0.0;
    return 2.0 * ml_atan(x / (1.0 + ml_sqrt(1.0 - x * x)));
}

static inline double ml_acos(double x) {
    if (x < -1.0 || x > 1.0) return 0.0 / 0.0;
    return (math_pi / 2.0) - ml_asin(x);
}

static inline double ml_acot(double x) {
    return (math_pi / 2.0) - ml_atan(x);
}

static inline double ml_atan2(double y, double x) {
    if (x == 0.0 && y == 0.0) return 0.0;
    if (x == 0.0) return (y > 0.0) ? math_pi / 2.0 : -math_pi / 2.0;
    double a = ml_atan(y / x);
    if (x < 0.0) return (y >= 0.0) ? a + math_pi : a - math_pi;
    return a;
}

#endif // MATHLIB_ML_TRIG_H
