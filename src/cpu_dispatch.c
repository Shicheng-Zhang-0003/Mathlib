#include "cpu_dispatch.h"
#include <stddef.h>

#if defined(__x86_64__) || defined(__i386__)
#include <immintrin.h>
#endif

/* MATHLIB_CLOSURE_P1_SIZE_T_INDEXING */

static void ml_matmul_scalar(const double* ML_RESTRICT A, const double* ML_RESTRICT B, double* ML_RESTRICT C, int N) {
    size_t n = (size_t)N;

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            double sum = 0.0;
            for (size_t k = 0; k < n; k++) {
                sum += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = sum;
        }
    }
}

#if ML_COMPILE_TIME_AVX2
ML_TARGET_AVX2
static void ml_matmul_avx2(const double* ML_RESTRICT A, const double* ML_RESTRICT B, double* ML_RESTRICT C, int N) {
    size_t n = (size_t)N;

    for (size_t i = 0; i < n; i++) {
        size_t j = 0;

        for (; j + 4 <= n; j += 4) {
            __m256d c_vec = _mm256_setzero_pd();

            for (size_t k = 0; k < n; k++) {
                __m256d a_vec = _mm256_broadcast_sd(&A[i * n + k]);
                __m256d b_vec = _mm256_loadu_pd(&B[k * n + j]);
                c_vec = _mm256_fmadd_pd(a_vec, b_vec, c_vec);
            }

            _mm256_storeu_pd(&C[i * n + j], c_vec);
        }

        for (; j < n; j++) {
            double sum = 0.0;
            for (size_t k = 0; k < n; k++) {
                sum += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = sum;
        }
    }
}
#endif

ML_API void ml_matmul(const double* ML_RESTRICT A, const double* ML_RESTRICT B, double* ML_RESTRICT C, int N) {
    if (ML_UNLIKELY(A == NULL || B == NULL || C == NULL || N <= 0)) return;

    size_t n = (size_t)N;

    /* Defensive guard against index-space overflow. */
    if (ML_UNLIKELY(n > ((size_t)-1) / n)) return;

#if ML_COMPILE_TIME_AVX2
    ml_matmul_avx2(A, B, C, N);
#else
    ml_matmul_scalar(A, B, C, N);
#endif
}

ML_API int ml_cpu_has_avx2(void) {
    return ML_COMPILE_TIME_AVX2;
}

ML_API int ml_cpu_has_fma(void) {
#if defined(__FMA__) && (defined(__x86_64__) || defined(__i386__))
    return 1;
#else
    return 0;
#endif
}

ML_API int ml_cpu_has_sse41(void) {
#if defined(__SSE4_1__) && (defined(__x86_64__) || defined(__i386__))
    return 1;
#else
    return 0;
#endif
}

ML_API int ml_cpu_has_neon(void) {
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
    return 1;
#else
    return 0;
#endif
}
