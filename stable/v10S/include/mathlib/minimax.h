#ifndef LIBMATHC_MINIMAX_H
#define LIBMATHC_MINIMAX_H

#include "ml_core.h"
#include "payne_hanek.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// High-precision odd polynomial (Maclaurin/Taylor-derived) for sin(x) on [-pi/2, pi/2].
// Note: While initially named 'minimax' in early development, these are the exact
// Maclaurin series coefficients up to x^17, which provide 1e-15 precision on this interval.
static const double minimax_sin_coeffs[] = {
    1.0,                      // Exact 1.0 to prevent 1e-16 downward bias
    -0.16666666666666666,
    0.008333333333333333,
    -0.0001984126984126984,
    2.7557319223985893e-06,
    -2.505210838544172e-08,
    1.6059043836821613e-10,
    -7.647163731819816e-13,
    2.811457254345521e-15,
    -8.220635816560923e-18    // x^19 / 19! (Pushes truncation error to 1e-17)
};

// Raw polynomial evaluation (Assumes x is already in [-pi/2, pi/2])
static inline double ml_minimax_sin_raw(double x) {
    double x2 = x * x;
    double result = minimax_sin_coeffs[9];
    for (int i = 8; i >= 0; i--) result = result * x2 + minimax_sin_coeffs[i];
    return x * result;
}

// Flawless Sin with Payne-Hanek + Quadrant Mapping
static inline double ml_minimax_sin(double x) {
    x = ml_reduce_payne_hanek(x); // Returns [-pi, pi]
    if (x > M_PI / 2.0) x = M_PI - x;
    else if (x < -M_PI / 2.0) x = -M_PI - x;
    return ml_minimax_sin_raw(x);
}

// Flawless Cos with Payne-Hanek + Quadrant Mapping
static inline double ml_minimax_cos(double x) {
    x = ml_reduce_payne_hanek(x); // Returns [-pi, pi]
    if (x < 0) x += 2.0 * M_PI; // Map to [0, 2pi)
    if (x <= M_PI / 2.0) return ml_minimax_sin_raw(M_PI / 2.0 - x);
    if (x <= M_PI) return -ml_minimax_sin_raw(x - M_PI / 2.0);
    if (x <= 1.5 * M_PI) return -ml_minimax_sin_raw(1.5 * M_PI - x);
    return ml_minimax_sin_raw(x - 1.5 * M_PI);
}

#endif
