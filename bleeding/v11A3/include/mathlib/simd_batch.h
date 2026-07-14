#ifndef MATHLIB_SIMD_BATCH_H
#define MATHLIB_SIMD_BATCH_H

#include "profiles.h"

#if defined(__AVX2__)
#include <immintrin.h>
typedef double ml_vec4d __attribute__((aligned(32), vector_size(32)));

static inline void ml_simd_batch_poly(const double* in, double* out) {
    __m256d v = _mm256_loadu_pd(in);
    __m256d res = _mm256_add_pd(_mm256_mul_pd(v, v), v);
    _mm256_storeu_pd(out, res);
}

static inline void ml_simd_batch_rsqrt(const double* in, double* out) {
    __m256d v = _mm256_loadu_pd(in);
    __m256i vi = _mm256_castpd_si256(v);
    __m256i magic = _mm256_set1_epi64x(0x5fe6ec85e7de30daLL);
    __m256i vi_half = _mm256_srli_epi64(vi, 1);
    __m256i yi = _mm256_sub_epi64(magic, vi_half);
    __m256d y = _mm256_castsi256_pd(yi);
    __m256d half = _mm256_set1_pd(0.5);
    __m256d three_half = _mm256_set1_pd(1.5);
    __m256d x_half = _mm256_mul_pd(v, half);
    __m256d y2 = _mm256_mul_pd(y, y);
    __m256d sub = _mm256_sub_pd(three_half, _mm256_mul_pd(x_half, y2));
    y = _mm256_mul_pd(y, sub);
    y2 = _mm256_mul_pd(y, y);
    sub = _mm256_sub_pd(three_half, _mm256_mul_pd(x_half, y2));
    y = _mm256_mul_pd(y, sub);
    _mm256_storeu_pd(out, y);
}

#else
// Scalar Fallback
static inline void ml_simd_batch_poly(const double* in, double* out) {
    for(int i=0; i<4; i++) out[i] = (in[i] * in[i]) + in[i];
}
static inline void ml_simd_batch_rsqrt(const double* in, double* out) {
    for(int i=0; i<4; i++) out[i] = (in[i] > 0.0) ? (1.0 / __builtin_sqrt(in[i])) : 0.0;
}
#endif

#endif
