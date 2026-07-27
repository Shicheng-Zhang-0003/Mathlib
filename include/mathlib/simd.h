#ifndef LIBMATHC_SIMD_H
#define LIBMATHC_SIMD_H

#include <stdint.h>
#include "ml_core.h"

#if defined(__AVX__) || defined(__AVX2__)
// GCC vector extension for 4-wide double precision (AVX)
typedef double ml_vec4 __attribute__((vector_size(32)));
#define ML_VEC4_AT(v, i) ((v)[i])
#define ML_VEC4_INIT(a,b,c,d) {a, b, c, d}

static inline ml_vec4 ml_vec4_add(ml_vec4 a, ml_vec4 b) { return a + b; }
static inline ml_vec4 ml_vec4_sub(ml_vec4 a, ml_vec4 b) { return a - b; }
static inline ml_vec4 ml_vec4_mul(ml_vec4 a, ml_vec4 b) { return a * b; }
static inline ml_vec4 ml_vec4_scale(ml_vec4 a, double s) { return a * (ml_vec4){s, s, s, s}; }

static inline double ml_vec4_dot(ml_vec4 a, ml_vec4 b) {
    ml_vec4 prod = a * b;
    return prod[0] + prod[1] + prod[2] + prod[3];
}

static inline double ml_vec4_mag(ml_vec4 a) {
    return ml_sqrt(ml_vec4_dot(a, a));
}

// --- Raw AVX Intrinsics (Path 3) ---
#include <immintrin.h>
static inline double ml_vec4_dot_avx(ml_vec4 a, ml_vec4 b) {
    __m256d va = _mm256_loadu_pd((double*)&a);
    __m256d vb = _mm256_loadu_pd((double*)&b);
    __m256d prod = _mm256_mul_pd(va, vb);
    __m128d vlow = _mm256_castpd256_pd128(prod);
    __m128d vhigh = _mm256_extractf128_pd(prod, 1);
    vlow = _mm_add_pd(vlow, vhigh);
    __m128d high64 = _mm_unpackhi_pd(vlow, vlow);
    return _mm_cvtsd_f64(_mm_add_sd(vlow, high64));
}

#else
// Scalar Fallback for ARM, RISC-V, and non-AVX x86
typedef struct { double d[4]; } ml_vec4;
#define ML_VEC4_AT(v, i) ((v).d[i])
#define ML_VEC4_INIT(a,b,c,d) {{a, b, c, d}}

static inline ml_vec4 ml_vec4_add(ml_vec4 a, ml_vec4 b) {
    return (ml_vec4){{a.d[0]+b.d[0], a.d[1]+b.d[1], a.d[2]+b.d[2], a.d[3]+b.d[3]}};
}
static inline ml_vec4 ml_vec4_sub(ml_vec4 a, ml_vec4 b) {
    return (ml_vec4){{a.d[0]-b.d[0], a.d[1]-b.d[1], a.d[2]-b.d[2], a.d[3]-b.d[3]}};
}
static inline ml_vec4 ml_vec4_mul(ml_vec4 a, ml_vec4 b) {
    return (ml_vec4){{a.d[0]*b.d[0], a.d[1]*b.d[1], a.d[2]*b.d[2], a.d[3]*b.d[3]}};
}
static inline ml_vec4 ml_vec4_scale(ml_vec4 a, double s) {
    return (ml_vec4){{a.d[0]*s, a.d[1]*s, a.d[2]*s, a.d[3]*s}};
}
static inline double ml_vec4_dot(ml_vec4 a, ml_vec4 b) {
    return a.d[0]*b.d[0] + a.d[1]*b.d[1] + a.d[2]*b.d[2] + a.d[3]*b.d[3];
}
static inline double ml_vec4_mag(ml_vec4 a) {
    return ml_sqrt(ml_vec4_dot(a, a));
}
#endif

#endif
