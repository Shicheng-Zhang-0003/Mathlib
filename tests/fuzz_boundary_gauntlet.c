#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include <time.h>
#include "ml_core.h"
#include "bitwise_fp.h"
#include "ml_trig.h"
#include "ml_exp_log.h"
#include "ml_fixed_point.h"
#include "ml_tensor.h"
#include "ml_linalg.h"
#include "ml_types.h"
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
    CHECK_NEAR(ml_sin(0.0), 0.0, 1e-15, "sin(0)");
    CHECK_NEAR(ml_sin(ML_PI/2.0), 1.0, 1e-15, "sin(pi/2)");
    CHECK_NEAR(ml_sin(ML_PI), 0.0, 1e-14, "sin(pi)");
    CHECK_NEAR(ml_sin(3.0*ML_PI/2.0), -1.0, 1e-14, "sin(3pi/2)");
    CHECK_NEAR(ml_sin(2.0*ML_PI), 0.0, 1e-14, "sin(2pi)");
    CHECK_NEAR(ml_cos(0.0), 1.0, 1e-15, "cos(0)");
    CHECK_NEAR(ml_cos(ML_PI/2.0), 0.0, 1e-15, "cos(pi/2)");
    CHECK_NEAR(ml_cos(ML_PI), -1.0, 1e-14, "cos(pi)");

    // Boundary crossing: 1 ULP away from pi/2
    double eps = 2.220446049250313e-16;
    CHECK_NEAR(ml_sin(ML_PI/2.0 - eps), ml_cos(eps), 1e-12, "sin(pi/2 - eps) == cos(eps)");
    CHECK_NEAR(ml_sin(ML_PI/2.0 + eps), ml_cos(eps), 1e-12, "sin(pi/2 + eps) == cos(eps)");
    CHECK_NEAR(ml_sin(-ML_PI/2.0 - eps), -ml_cos(eps), 1e-12, "sin(-pi/2 - eps) == -cos(eps)");
}

void test_exp_log_boundaries() {
    printf("--- Exp/Log Exact Boundaries ---\n");
    CHECK_NEAR(ml_log(1.0), 0.0, 1e-15, "log(1)");
    CHECK_NEAR(ml_log(ML_E), 1.0, 1e-14, "log(e)");
    CHECK_NEAR(ml_exp(0.0), 1.0, 1e-15, "exp(0)");
    CHECK_NEAR(ml_exp(1.0), ML_E, 1e-14, "exp(1)");
    CHECK(ml_isinf(ml_log(0.0)) && ml_log(0.0) < 0, "log(0) == -inf");
    CHECK(ml_isinf(ml_exp(1000.0)) && ml_exp(1000.0) > 0, "exp(1000) == +inf");
    CHECK(ml_exp(-1000.0) < 5e-324, "exp(-1000) underflows to zero or subnormal");
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
    char scratchpad[4096]; // FIXED: 1MB stack bomb reduced to 4KB
    ml_workspace_t ws;
    ml_workspace_init(&ws, scratchpad, sizeof(scratchpad));
    int n = 4;
    double A_data[16];
    double b_data[4] = {1.0, 1.0, 1.0, 1.0};
    double x_data[4] = {0};
    double residual[4] = {0};

    for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
            A_data[i*n+j] = 1.0 / (i + j + 1.0);
        }
    }
    ml_tensor_view_t A_view = ml_tensor_view(A_data, n, n);
    int status = ml_solve(A_view, b_data, x_data, &ws);
    CHECK(status == ML_SUCCESS, "Hilbert solve status");

    // FIXED: Actually validate the solution using residual ||Ax - b||
    ml_matvec(A_view, x_data, residual);
    double max_err = 0.0;
    for(int i=0; i<n; i++) {
        double err = ml_fabs(residual[i] - b_data[i]);
        if (err > max_err) max_err = err;
    }
    // Hilbert n=4 condition number is ~15,000. 1e-10 is a strict but fair bound.
    CHECK(max_err < 1e-10, "Hilbert residual ||Ax - b|| < 1e-10");
}

void test_fft_parseval() {
    printf("--- FFT Parseval's Theorem (Energy Conservation) ---\n");
    srand(42); // FIXED: Deterministic seed for exact reproduction
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
    ml_fft_execute(freq_domain, n);
    double freq_energy = 0.0;
    for(int i=0; i<n; i++) {
        freq_energy += (freq_domain[i].real * freq_domain[i].real + freq_domain[i].imag * freq_domain[i].imag);
    }
    freq_energy /= n;
    // FIXED: Tightened from 1e-6 (useless) to 1e-14 (physical limit of double precision)
    CHECK(ml_fabs(time_energy - freq_energy) < time_energy * 1e-13, "Parseval's Theorem (Time Energy == Freq Energy)");
}

int main() {
    printf("=========================================================\n");
    printf("   MATHLIB v11S: BOUNDARY & INVARIANT GAUNTLET\n");
    printf("=========================================================\n");
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
