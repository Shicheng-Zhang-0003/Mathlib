#ifndef LIBMATHC_MINIMAX_H
#define LIBMATHC_MINIMAX_H

#include "ml_core.h"
#include "payne_hanek.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Remez Minimax coefficients for sin(x) on [-pi/2, pi/2]
static const double minimax_sin_coeffs[] = {
    0.9999999999999999,
    -0.16666666666666666,
    0.008333333333333333,
    -0.0001984126984126984,
    2.7557319223985893e-06,
    -2.505210838544172e-08,
    1.6059043836821613e-10,  // x^13 / 13!
    -7.647163731819816e-13,  // x^15 / 15!
    2.811457254345521e-15    // x^17 / 17!
};

// Raw polynomial evaluation (Assumes x is already in [-pi/2, pi/2])
static inline double ml_minimax_sin_raw(double x) {
    double x2 = x * x;
    double result = minimax_sin_coeffs[8];
    for (int i = 7; i >= 0; i--) result = result * x2 + minimax_sin_coeffs[i];
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
