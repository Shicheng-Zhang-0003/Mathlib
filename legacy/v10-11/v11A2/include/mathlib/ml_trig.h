#ifndef MATHLIB_ML_TRIG_H
#define MATHLIB_ML_TRIG_H

#include "ml_core.h"
#include "internal/minimax.h"
#include "internal/cordic.h"
#include "fixed_point.h" // Required for ml_q16_16_t and ml_cordic_sincos_fixed

#ifndef math_pi
#define math_pi 3.14159265358979323846
#endif

// --- Polymorphic API (Simultaneous Profile Access) ---
static inline double ml_sin_scientific(double x) { return ml_minimax_sin(x); }
static inline double ml_cos_scientific(double x) { return ml_minimax_cos(x); }
static inline double ml_sin_graphics(double x) { return ml_minimax_sin(x); }
static inline double ml_cos_graphics(double x) { return ml_minimax_cos(x); }

static inline double ml_sin_fixed(double x) {
    x = ml_fmod(x, 2.0 * 3.14159265358979323846);
    ml_q16_16_t f_in = (ml_q16_16_t)(x * 65536.0);
    ml_q16_16_t s, c;
    ml_cordic_sincos_fixed(f_in, &s, &c);
    return (double)s / 65536.0;
}

static inline double ml_cos_fixed(double x) {
    x = ml_fmod(x, 2.0 * 3.14159265358979323846);
    ml_q16_16_t f_in = (ml_q16_16_t)(x * 65536.0);
    ml_q16_16_t s, c;
    ml_cordic_sincos_fixed(f_in, &s, &c);
    return (double)c / 65536.0;
}

// Profile-Routed Trigonometry (The v11A2 Way)
static inline double ml_sin(double x) {
    if (ml_isnan(x) || ml_isinf(x)) return 0.0/0.0;
#if defined(MATHLIB_PROFILE_EMBEDDED)
    return ml_sin_fixed(x);
#else
    return ml_minimax_sin(x);
#endif
}

static inline double ml_cos(double x) {
    if (ml_isnan(x) || ml_isinf(x)) return 0.0/0.0;
#if defined(MATHLIB_PROFILE_EMBEDDED)
    return ml_cos_fixed(x);
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

static inline double ml_acot(double x) { return (math_pi / 2.0) - ml_atan(x); }

static inline double ml_atan2(double y, double x) {
    if (x == 0.0 && y == 0.0) return 0.0;
    if (x == 0.0) return (y > 0.0) ? math_pi / 2.0 : -math_pi / 2.0;
    double a = ml_atan(y / x);
    if (x < 0.0) return (y >= 0.0) ? a + math_pi : a - math_pi;
    return a;
}

#endif // MATHLIB_ML_TRIG_H
