#ifndef LIBMATHC_SIMD_H
#define LIBMATHC_SIMD_H

// GCC vector extension for 4-wide double precision (AVX)
typedef double ml_vec4 __attribute__((vector_size(32)));

inline ml_vec4 ml_vec4_add(ml_vec4 a, ml_vec4 b) { return a + b; }
inline ml_vec4 ml_vec4_sub(ml_vec4 a, ml_vec4 b) { return a - b; }
inline ml_vec4 ml_vec4_mul(ml_vec4 a, ml_vec4 b) { return a * b; }
inline ml_vec4 ml_vec4_scale(ml_vec4 a, double s) { return a * (ml_vec4){s, s, s, s}; }

inline double ml_vec4_dot(ml_vec4 a, ml_vec4 b) {
    ml_vec4 prod = a * b;
    return prod[0] + prod[1] + prod[2] + prod[3];
}

inline double ml_vec4_mag(ml_vec4 a) {
    double dot = ml_vec4_dot(a, a);
    double res;
    __asm__ ("sqrtsd %1, %0" : "=x" (res) : "x" (dot));
    return res;
}


// --- Raw AVX Intrinsics (Path 3) ---
#include <immintrin.h>

inline double ml_vec4_dot_avx(ml_vec4 a, ml_vec4 b) {
    __m256d va = _mm256_loadu_pd((double*)&a);
    __m256d vb = _mm256_loadu_pd((double*)&b);
    __m256d prod = _mm256_mul_pd(va, vb);

    // Horizontal add across the 256-bit register
    __m128d vlow = _mm256_castpd256_pd128(prod);
    __m128d vhigh = _mm256_extractf128_pd(prod, 1);
    vlow = _mm_add_pd(vlow, vhigh);
    __m128d high64 = _mm_unpackhi_pd(vlow, vlow);
    return _mm_cvtsd_f64(_mm_add_sd(vlow, high64));
}

#endif
