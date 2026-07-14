#ifndef MATHLIB_CPU_DISPATCH_H
#define MATHLIB_CPU_DISPATCH_H

typedef struct {
    int avx2_supported;
    int fma_supported;
    int sse41_supported;
    int neon_supported;
    int initialized;
} ml_cpu_caps_t;

extern ml_cpu_caps_t g_ml_cpu_caps;

void ml_init_cpu_dispatch(void);

typedef void (*ml_matmul_func_t)(const double* A, const double* B, double* C, int N);
extern ml_matmul_func_t g_ml_matmul;

void ml_matmul_scalar(const double* A, const double* B, double* C, int N);
void ml_matmul_avx2(const double* A, const double* B, double* C, int N);

static inline void ml_matmul(const double* A, const double* B, double* C, int N) {
    if (!g_ml_cpu_caps.initialized) ml_init_cpu_dispatch();
    g_ml_matmul(A, B, C, N);
}

#endif
