#ifndef LIBMATHC_TRIGONOMETRY_H
#define LIBMATHC_TRIGONOMETRY_H

#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define math_pi M_PI
inline double arctangent(double x);

// --- Abstraction Layer for Range Reduction ---
// Future-proofed: Swap the inside of this function with IEEE 754 bit-masking later
inline double reduce_angle(double x) {
    x = fmod(x, 2.0 * math_pi);
    if (x > math_pi) x -= 2.0 * math_pi;
    if (x < -math_pi) x += 2.0 * math_pi;
    return x;
}

inline double sine(double x) {
    x = reduce_angle(x);
    double result = x;
    double term = x;
    double x2 = x * x;
    for (int step = 3; step <= 21; step += 2) {
        term *= -x2 / ((step - 1) * step);
        result += term;
    }
    return result;
}

inline double cosine(double x) {
    x = reduce_angle(x);
    double result = 1.0;
    double term = 1.0;
    double x2 = x * x;
    for (int step = 2; step <= 20; step += 2) {
        term *= -x2 / ((step - 1) * step);
        result += term;
    }
    return result;
}

inline double tangent(double x) {
    double c = cosine(x);
    if (c == 0.0) return 0.0 / 0.0; // NaN for asymptotes
    return sine(x) / c;
}

inline double cosecant(double x) { double s = sine(x); return s == 0.0 ? 0.0/0.0 : 1.0 / s; }
inline double secant(double x) { double c = cosine(x); return c == 0.0 ? 0.0/0.0 : 1.0 / c; }
inline double cotangent(double x) { double s = sine(x); return s == 0.0 ? 0.0/0.0 : cosine(x) / s; }

inline double arcsine(double x) {
    if (x < -1.0 || x > 1.0) return 0.0 / 0.0;
    return arctangent(x / sqrt(1.0 - x * x));
}

inline double arccosine(double x) {
    if (x < -1.0 || x > 1.0) return 0.0 / 0.0;
    return (math_pi / 2.0) - arcsine(x);
}

inline double arctangent(double x) {
    if (x > 1.0) return (math_pi / 2.0) - arctangent(1.0 / x);
    if (x < -1.0) return -(math_pi / 2.0) - arctangent(1.0 / x);
    if (x > 0.5) return (math_pi / 4.0) + arctangent((x - 1.0) / (x + 1.0));
    if (x < -0.5) return -(math_pi / 4.0) + arctangent((x + 1.0) / (1.0 - x));

    double result = x;
    double term = x;
    double x2 = x * x;
    for (int step = 3; step <= 21; step += 2) {
        term *= -x2;
        result += term / step;
    }
    return result;
}

inline double arccosecant(double x) { return (x <= -1.0 || x >= 1.0) ? arcsine(1.0 / x) : 0.0/0.0; }
inline double arcsecant(double x) { return (x <= -1.0 || x >= 1.0) ? arccosine(1.0 / x) : 0.0/0.0; }
inline double arccotangent(double x) {
    if (x == 0.0) return math_pi / 2.0;
    return (x > 0.0) ? arctangent(1.0 / x) : math_pi + arctangent(1.0 / x);
}

#endif
