#ifndef MATHLIB_SIMD_BATCH_H
#define MATHLIB_SIMD_BATCH_H

#include "profiles.h"
#include "ml_core.h"
#include "bitwise_fp.h"
#include "fast_math.h"

/* MATHLIB_CLOSURE_P1_SIMD_BATCH_GUARD */
/* MATHLIB_REDP2_P1_1_SIMD_BATCH_SCALAR_FALLBACK */

#if defined(__AVX2__)
#include <immintrin.h>

typedef double ml_vec4d __attribute__((aligned(32), vector_size(32)));

static inline void ml_simd_batch_poly(const double* in, double* out) {
    __m256d v = _mm256_loadu_pd(in);
    __m256d res = _mm256_add_pd(_mm256_mul_pd(v, v), v);
    _mm256_storeu_pd(out, res);
}

static inline void ml_simd_batch_rsqrt(const double* in, double* out) {
    /*
     * MATHLIB_CLOSURE_P1:
     * Exceptional inputs fall back to scalar semantics for the whole batch.
     * The AVX2 fast path is used only for positive finite normal inputs.
     *
     * MATHLIB_REDP2_P1_1:
     * The scalar fallback now uses ml_fast_rsqrt() per lane, so batch
     * exceptional behavior matches the scalar fast-math contract exactly:
     *
     *   NaN             -> NaN
     *   negative finite -> NaN
     *   +inf            -> +0
     *   +0              -> +inf
     *   -0              -> -inf
     *   positive finite -> approximate rsqrt
     */
    for (int i = 0; i < 4; i++) {
        if (ml_isnan(in[i]) || ml_isinf(in[i]) ||
            ml_is_subnormal(in[i]) || !(in[i] > 0.0)) {
            for (int j = 0; j < 4; j++) {
                out[j] = ml_fast_rsqrt(in[j]);
            }
            return;
        }
    }

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

/* Scalar Fallback */

static inline void ml_simd_batch_poly(const double* in, double* out) {
    for (int i = 0; i < 4; i++) {
        out[i] = (in[i] * in[i]) + in[i];
    }
}

static inline void ml_simd_batch_rsqrt(const double* in, double* out) {
    /*
     * MATHLIB_REDP2_P1_1:
     * Use the scalar fast-math contract directly.
     */
    for (int i = 0; i < 4; i++) {
        out[i] = ml_fast_rsqrt(in[i]);
    }
}

#endif

#endif /* MATHLIB_SIMD_BATCH_H */
