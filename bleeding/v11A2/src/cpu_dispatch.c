#include "cpu_dispatch.h"

#if defined(__x86_64__) || defined(__i386__)
#include <immintrin.h>
#endif

#if defined(__aarch64__)
#include <sys/auxv.h>
#ifndef HWCAP_ASIMD
#define HWCAP_ASIMD (1 << 1)
#endif
#ifndef HWCAP_FPHP
#define HWCAP_FPHP (1 << 9)
#endif
#endif

ml_cpu_caps_t g_ml_cpu_caps = {0};
ml_matmul_func_t g_ml_matmul = ml_matmul_scalar;

void ml_matmul_scalar(const double* A, const double* B, double* C, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double sum = 0.0;
            for (int k = 0; k < N; k++) {
                sum += A[i * N + k] * B[k * N + j];
            }
            C[i * N + j] = sum;
        }
    }
}

#if defined(__x86_64__) || defined(__i386__)
__attribute__((target("avx2,fma")))
void ml_matmul_avx2(const double* A, const double* B, double* C, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j += 4) {
            __m256d c_vec = _mm256_setzero_pd();
            for (int k = 0; k < N; k++) {
                __m256d a_vec = _mm256_broadcast_sd(&A[i * N + k]);
                __m256d b_vec = _mm256_loadu_pd(&B[k * N + j]);
                c_vec = _mm256_fmadd_pd(a_vec, b_vec, c_vec);
            }
            _mm256_storeu_pd(&C[i * N + j], c_vec);
        }
    }
}
#else
void ml_matmul_avx2(const double* A, const double* B, double* C, int N) {
    ml_matmul_scalar(A, B, C, N);
}
#endif

void ml_init_cpu_dispatch(void) {
    if (g_ml_cpu_caps.initialized) return;

#if defined(__x86_64__) || defined(__i386__)
    g_ml_cpu_caps.avx2_supported = __builtin_cpu_supports("avx2");
    g_ml_cpu_caps.fma_supported = __builtin_cpu_supports("fma");
    g_ml_cpu_caps.sse41_supported = __builtin_cpu_supports("sse4.1");
    g_ml_cpu_caps.neon_supported = 0;
#elif defined(__aarch64__)
    long hwcaps = getauxval(AT_HWCAP);
    g_ml_cpu_caps.neon_supported = (hwcaps & HWCAP_ASIMD) ? 1 : 0;
    g_ml_cpu_caps.avx2_supported = 0;
    g_ml_cpu_caps.fma_supported = (hwcaps & HWCAP_FPHP) ? 1 : 0;
#else
    g_ml_cpu_caps.avx2_supported = 0;
    g_ml_cpu_caps.fma_supported = 0;
    g_ml_cpu_caps.sse41_supported = 0;
    g_ml_cpu_caps.neon_supported = 0;
#endif

    if (g_ml_cpu_caps.avx2_supported && g_ml_cpu_caps.fma_supported) {
        g_ml_matmul = ml_matmul_avx2;
    } else {
        g_ml_matmul = ml_matmul_scalar;
    }

    g_ml_cpu_caps.initialized = 1;
}

__attribute__((constructor))
static void ml_auto_init_cpu_dispatch(void) {
    ml_init_cpu_dispatch();
}
