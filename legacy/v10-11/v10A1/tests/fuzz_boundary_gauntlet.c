
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include <string.h>
#include <time.h>

#include "ml_core.h"
#include "bitwise_fp.h"
#include "trigonometry.h"
#include "exponential.h"
#include "fixed_point.h"
#include "tensor.h"
#include "linalg_v10.h"
#include "fft.h"
#include "ml_complex.h"
#include "fast_math.h"

int passed = 0;
int failed = 0;

#define CHECK(cond, msg) do { \
    if (cond) { passed++; } \
    else { failed++; printf("FAIL: %s (Line %d)\n", msg, __LINE__); } \
} while(0)

#define CHECK_NEAR(a, b, eps, msg) CHECK(ml_fabs((double)(a) - (double)(b)) < (eps), msg)

void test_trig_boundaries() {
    printf("--- Trig Exact Boundaries ---\n");
    double pi = 3.14159265358979323846;

    CHECK_NEAR(sine(0.0), 0.0, 1e-15, "sin(0)");
    CHECK_NEAR(sine(pi/2.0), 1.0, 1e-15, "sin(pi/2)");
    CHECK_NEAR(sine(pi), 0.0, 1e-14, "sin(pi)");
    CHECK_NEAR(sine(3.0*pi/2.0), -1.0, 1e-14, "sin(3pi/2)");
    CHECK_NEAR(sine(2.0*pi), 0.0, 1e-14, "sin(2pi)");

    CHECK_NEAR(cosine(0.0), 1.0, 1e-15, "cos(0)");
    CHECK_NEAR(cosine(pi/2.0), 0.0, 1e-15, "cos(pi/2)");
    CHECK_NEAR(cosine(pi), -1.0, 1e-14, "cos(pi)");

    // Boundary crossing: just inside and outside the [-pi/2, pi/2] minimax domain
    double eps = 1e-9;
    CHECK_NEAR(sine(pi/2.0 - eps), cosine(eps), 1e-12, "sin(pi/2 - eps) == cos(eps)");
    CHECK_NEAR(sine(pi/2.0 + eps), cosine(eps), 1e-12, "sin(pi/2 + eps) == cos(eps)");
    CHECK_NEAR(sine(-pi/2.0 - eps), -cosine(eps), 1e-12, "sin(-pi/2 - eps) == -cos(eps)");
}

void test_exp_log_boundaries() {
    printf("--- Exp/Log Exact Boundaries ---\n");
    double e = 2.71828182845904523536;

    CHECK_NEAR(logarithm(1.0), 0.0, 1e-15, "log(1)");
    CHECK_NEAR(logarithm(e), 1.0, 1e-14, "log(e)");
    CHECK_NEAR(exponential(0.0), 1.0, 1e-15, "exp(0)");
    CHECK_NEAR(exponential(1.0), e, 1e-14, "exp(1)");

    CHECK(ml_isinf(logarithm(0.0)) && logarithm(0.0) < 0, "log(0) == -inf");
    CHECK(exponential(1000.0) == (1.0/0.0), "exp(1000) == inf");
    CHECK(exponential(-1000.0) == 0.0, "exp(-1000) == 0");
}

void test_fixed_point_saturation() {
    printf("--- Fixed-Point Saturation ---\n");

    ml_q16_16_t max_val = 2147483647;
    ml_q16_16_t min_val = (-2147483647 - 1);

    CHECK(ml_fixed_mul(max_val, max_val) == 2147483647, "MAX * MAX saturates to MAX");
    CHECK(ml_fixed_mul(min_val, max_val) == (-2147483647 - 1), "MIN * MAX saturates to MIN");
    CHECK(ml_fixed_mul(min_val, min_val) == 2147483647, "MIN * MIN saturates to MAX");

    ml_q16_16_t two = 2 << 16;
    ml_q16_16_t three = 3 << 16;
    CHECK(ml_fixed_mul(two, three) == (6 << 16), "2.0 * 3.0 == 6.0");
}

void test_tensor_hilbert() {
    printf("--- Tensor Ill-Conditioned (Hilbert) ---\n");
    char scratchpad[1024 * 1024];
    ml_workspace_t ws = { scratchpad, sizeof(scratchpad), 0 };

    int n = 4;
    double A_data[16];
    double b_data[4] = {1.0, 1.0, 1.0, 1.0};
    double x_data[4] = {0};

    // 4x4 Hilbert Matrix (Notoriously ill-conditioned)
    for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
            A_data[i*n+j] = 1.0 / (i + j + 1.0);
        }
    }

    ml_tensor_view_t A_view = ml_tensor_view(A_data, n, n);
    int status = ml_solve_v10(A_view, b_data, x_data, &ws);

    CHECK(status == 0, "Hilbert solve status");
    int has_nan = 0;
    for(int i=0; i<n; i++) if (ml_isnan(x_data[i]) || ml_isinf(x_data[i])) has_nan = 1;
    CHECK(!has_nan, "Hilbert solution has no NaN/Inf");
}

void test_fft_parseval() {
    printf("--- FFT Parseval's Theorem (Energy Conservation) ---\n");
    int n = 64;
    cplx time_domain[64];
    cplx freq_domain[64];

    double time_energy = 0.0;
    for(int i=0; i<n; i++) {
        double val = (double)rand() / RAND_MAX * 10.0 - 5.0;
        time_domain[i] = (cplx){val, 0.0};
        time_energy += val * val;
    }

    for(int i=0; i<n; i++) freq_domain[i] = time_domain[i];
    fft_execute(freq_domain, n);

    double freq_energy = 0.0;
    for(int i=0; i<n; i++) {
        freq_energy += (freq_domain[i].real * freq_domain[i].real + freq_domain[i].imag * freq_domain[i].imag);
    }
    freq_energy /= n;

    CHECK_NEAR(time_energy, freq_energy, 1e-6, "Parseval's Theorem (Time Energy == Freq Energy)");
}

int main() {
    srand(time(NULL));
    printf("=========================================================\n");
    printf("   MATHLIB v1.0: BOUNDARY & INVARIANT GAUNTLET\n");
    printf("=========================================================\n\n");

    test_trig_boundaries();
    test_exp_log_boundaries();
    test_fixed_point_saturation();
    test_tensor_hilbert();
    test_fft_parseval();

    printf("\n=========================================================\n");
    printf("BOUNDARY SUMMARY: %d passed, %d failed\n", passed, failed);
    printf("=========================================================\n");

    return failed > 0 ? 1 : 0;
}
