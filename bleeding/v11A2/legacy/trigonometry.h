#ifndef LIBMATHC_TRIGONOMETRY_H
#define LIBMATHC_TRIGONOMETRY_H

#include "ml_core.h"
#include "internal/minimax.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846



// Safely placed outside of any M_PI macro guards

#endif
#define math_pi M_PI

// Bypass the naive Taylor series and ml_fmod entirely.
// Route directly to the Remez Minimax polynomials which use
// the bulletproof Payne-Hanek zero-loss range reduction.

static inline double ml_atan(double x);

static inline double ml_sin(double x) {
    return ml_minimax_sin(x);
}

static inline double ml_cos(double x) {
    return ml_minimax_cos(x);
}

static inline double ml_tan(double x) {
    return ml_sin(x) / ml_cos(x);
}

static inline double ml_asin(double x) {
    if (x < -1.0 || x > 1.0) return 0.0 / 0.0;
    return 2.0 * ml_atan(x / (1.0 + ml_sqrt(1.0 - x * x))); // Half-angle identity prevents cancellation
}

static inline double ml_acos(double x) {
    if (x < -1.0 || x > 1.0) return 0.0 / 0.0;
    return (math_pi / 2.0) - ml_asin(x);
}

static inline double ml_atan(double x) {
    if (x > 1.0) return (math_pi / 2.0) - ml_atan(1.0 / x);
    if (x < -1.0) return -(math_pi / 2.0) - ml_atan(1.0 / x);
    if (x > 0.5) return (math_pi / 4.0) + ml_atan((x - 1.0) / (x + 1.0));
    if (x < -0.5) return -(math_pi / 4.0) + ml_atan((x + 1.0) / (1.0 - x));

    double result = x;
    double term = x;
    double x2 = x * x;
    for (int i = 3; i <= 21; i += 2) {
        term *= -x2;
        result += term / i;
    }
    return result;
}

static inline double ml_atan2(double y, double x) {
    if (x == 0.0 && y == 0.0) return 0.0;
    if (x == 0.0) return (y > 0.0) ? math_pi / 2.0 : -math_pi / 2.0;
    double a = ml_atan(y / x);
    if (x < 0.0) {
        return (y >= 0.0) ? a + math_pi : a - math_pi;
    }
    return a;
}



// Safely placed outside of any M_PI macro guards


// Safely placed outside of any M_PI macro guards
static inline double ml_acot(double x) {
    return (math_pi / 2.0) - ml_atan(x);
}

#endif
