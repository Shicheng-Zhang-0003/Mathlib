#ifndef MATHLIB_SIMD_BATCH_H
#define MATHLIB_SIMD_BATCH_H

#include "profiles.h"

typedef double ml_vec4d __attribute__((aligned(32), vector_size(32)));

static inline void ml_simd_batch_poly(const double* in, double* out) {
    ml_vec4d v = *(ml_vec4d*)in;
    ml_vec4d res = v * v + v;
    *(ml_vec4d*)out = res;
}

static inline void ml_simd_batch_rsqrt(const double* in, double* out) {
    ml_vec4d v = *(ml_vec4d*)in;
    ml_vec4d half = {0.5, 0.5, 0.5, 0.5};
    ml_vec4d three_half = {1.5, 1.5, 1.5, 1.5};
    ml_vec4d y = {1.0/in[0], 1.0/in[1], 1.0/in[2], 1.0/in[3]};
    y = y * (three_half - (v * half * y * y));
    *(ml_vec4d*)out = y;
}
#endif
