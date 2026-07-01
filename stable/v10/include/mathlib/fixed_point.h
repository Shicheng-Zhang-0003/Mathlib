#ifndef MATHLIB_V10_FIXED_POINT_H
#define MATHLIB_V10_FIXED_POINT_H

#include <stdint.h>

// Q16.16 Fixed-Point Format (16 bits integer, 16 bits fractional)
typedef int32_t ml_q16_16_t;

#define ML_FIXED_PI 205887 // 3.14159... * 65536
#define ML_FIXED_HALF_PI 102943
#define ML_FIXED_TWO_PI 411774

// Fixed-point multiplication
static inline ml_q16_16_t ml_fixed_mul(ml_q16_16_t a, ml_q16_16_t b) {
    int64_t res = ((int64_t)a * b) >> 16;
    // Saturation arithmetic to prevent wrap-around overflow in ML/DSP pipelines
    if (res > INT32_MAX) return INT32_MAX;
    if (res < INT32_MIN) return INT32_MIN;
    return (ml_q16_16_t)res;
}

// Precomputed CORDIC atan table in Q16.16 format
static const ml_q16_16_t fixed_cordic_atan[] = {
    51471, 30385, 16061, 8152, 4093, 2048, 1024, 512,
    256, 128, 64, 32, 16, 8, 4, 2, 1
};
#define ML_FIXED_CORDIC_GAIN 39797 // 0.60725... * 65536

// True Bare-Metal CORDIC (Zero floats, zero FPU dependencies)
static inline void ml_cordic_sincos_fixed(ml_q16_16_t theta, ml_q16_16_t *sin_out, ml_q16_16_t *cos_out) {
    // O(1) Integer Modulo range reduction to [-PI, PI]
    theta = theta % ML_FIXED_TWO_PI;
    if (theta > ML_FIXED_PI) theta -= ML_FIXED_TWO_PI;
    if (theta < -ML_FIXED_PI) theta += ML_FIXED_TWO_PI;

    int negate_cos = 0;
    if (theta > ML_FIXED_HALF_PI) {
        theta = ML_FIXED_PI - theta;
        negate_cos = 1;
    } else if (theta < -ML_FIXED_HALF_PI) {
        theta = -ML_FIXED_PI - theta;
        negate_cos = 1;
    }

    ml_q16_16_t x = ML_FIXED_CORDIC_GAIN;
    ml_q16_16_t y = 0;
    ml_q16_16_t z = theta;

    for (int i = 0; i < 16; i++) {
        ml_q16_16_t x_new, y_new;
        if (z >= 0) {
            x_new = x - (y >> i); // Pure bitwise shift instead of division!
            y_new = y + (x >> i);
            z -= fixed_cordic_atan[i];
        } else {
            x_new = x + (y >> i);
            y_new = y - (x >> i);
            z += fixed_cordic_atan[i];
        }
        x = x_new;
        y = y_new;
    }

    if (negate_cos) x = -x;
    *cos_out = x;
    *sin_out = y;
}
#endif
