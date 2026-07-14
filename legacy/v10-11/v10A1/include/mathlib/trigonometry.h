#ifndef LIBMATHC_TRIGONOMETRY_H
#define LIBMATHC_TRIGONOMETRY_H

#include "ml_core.h"
#include "minimax.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846



// Safely placed outside of any M_PI macro guards

#endif
#define math_pi M_PI

// Bypass the naive Taylor series and ml_fmod entirely.
// Route directly to the Remez Minimax polynomials which use
// the bulletproof Payne-Hanek zero-loss range reduction.

static inline double arctangent(double x);

static inline double sine(double x) {
    return ml_minimax_sin(x);
}

static inline double cosine(double x) {
    return ml_minimax_cos(x);
}

static inline double tangent(double x) {
    return sine(x) / cosine(x);
}

static inline double arcsine(double x) {
    if (x < -1.0 || x > 1.0) return 0.0 / 0.0;
    return 2.0 * arctangent(x / (1.0 + ml_sqrt(1.0 - x * x))); // Half-angle identity prevents cancellation
}

static inline double arccosine(double x) {
    if (x < -1.0 || x > 1.0) return 0.0 / 0.0;
    return (math_pi / 2.0) - arcsine(x);
}

static inline double arctangent(double x) {
    if (x > 1.0) return (math_pi / 2.0) - arctangent(1.0 / x);
    if (x < -1.0) return -(math_pi / 2.0) - arctangent(1.0 / x);
    if (x > 0.5) return (math_pi / 4.0) + arctangent((x - 1.0) / (x + 1.0));
    if (x < -0.5) return -(math_pi / 4.0) + arctangent((x + 1.0) / (1.0 - x));

    double result = x;
    double term = x;
    double x2 = x * x;
    for (int i = 3; i <= 21; i += 2) {
        term *= -x2;
        result += term / i;
    }
    return result;
}

static inline double arctangent2(double y, double x) {
    if (x == 0.0 && y == 0.0) return 0.0;
    if (x == 0.0) return (y > 0.0) ? math_pi / 2.0 : -math_pi / 2.0;
    double a = arctangent(y / x);
    if (x < 0.0) {
        return (y >= 0.0) ? a + math_pi : a - math_pi;
    }
    return a;
}



// Safely placed outside of any M_PI macro guards


// Safely placed outside of any M_PI macro guards
static inline double arccotangent(double x) {
    return (math_pi / 2.0) - arctangent(x);
}

#endif
