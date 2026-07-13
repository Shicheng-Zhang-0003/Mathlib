#ifndef MATHLIB_V10_SIMD_BARE_METAL_H
#define MATHLIB_V10_SIMD_BARE_METAL_H

#if defined(__AVX__) || defined(__AVX2__)
#include <immintrin.h>

// True Bare-Metal AVX2 Fast RSqrt
// Uses integer bit-hacks on the 256-bit vector register to generate the
// initial guess in 1 cycle, completely eliminating the scalar division bottleneck.
static inline __m256d ml_avx2_fast_rsqrt(__m256d v) {
    // 1. Cast double vector to 64-bit integer vector (Zero-cost cast)
    __m256i vi = _mm256_castpd_si256(v);
    // 2. The Magic Constant for 64-bit doubles, broadcasted to all 4 lanes
    __m256i magic = _mm256_set1_epi64x(0x5fe6ec85e7de30daLL);
    // 3. Integer shift and subtract (The Quake III hack, vectorized!)
    __m256i vi_half = _mm256_srli_epi64(vi, 1);
    __m256i yi = _mm256_sub_epi64(magic, vi_half);
    // 4. Cast back to double vector
    __m256d y = _mm256_castsi256_pd(yi);
    // 5. One iteration of Vectorized Newton-Raphson
    __m256d half = _mm256_set1_pd(0.5);
    __m256d three_half = _mm256_set1_pd(1.5);
    // y = y * (1.5 - (x * 0.5 * y * y))
    __m256d y2 = _mm256_mul_pd(y, y);
    __m256d x_half = _mm256_mul_pd(v, half);
    __m256d sub = _mm256_sub_pd(three_half, _mm256_mul_pd(x_half, y2));
    return _mm256_mul_pd(y, sub);
}

#endif // __AVX__ || __AVX2__

#endif // MATHLIB_V10_SIMD_BARE_METAL_H
