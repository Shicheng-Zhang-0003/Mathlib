
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

// Include the entire v10 engine
#include "ml_core.h"
#include "bitwise_fp.h"
#include "combinatorics.h"
#include "quadratics.h"
#include "integral.h"
#include "trigonometry.h"
#include "exponential.h"
#include "numerical.h"
#include "polynomial.h"
#include "complex.h"
#include "linear_algebra.h"
#include "statistics.h"
#include "ode.h"
#include "optimization.h"
#include "fft.h"
#include "fast_math.h"
#include "error_free.h"
#include "cordic.h"
#include "payne_hanek.h"
#include "minimax.h"
#include "matrix.h"
#include "ieee754.h"
#include "simd.h"
#include "quaternion.h"
#include "tensor.h"
#include "linalg_v10.h"
#include "fixed_point.h"
#include "simd_bare_metal.h"
#include "simd_batch.h"
#include "profiles.h"

int passed = 0;
int failed = 0;

#define CHECK(cond, msg) do { \
    if (cond) { passed++; } \
    else { failed++; printf("FAIL: %s (Line %d)\n", msg, __LINE__); } \
} while(0)

#define CHECK_NEAR(a, b, eps, msg) CHECK(fabs((double)(a) - (double)(b)) < (eps), msg)
#define CHECK_NAN(a, msg) CHECK(isnan(a), msg)
#define CHECK_INF(a, msg) CHECK(isinf(a), msg)

double rand_double() {
    int type = rand() % 5;
    if (type == 0) return (double)rand() / RAND_MAX * 2000.0 - 1000.0; // Normal
    if (type == 1) return (double)rand() / RAND_MAX * 1e-15; // Tiny
    if (type == 2) return (double)rand() / RAND_MAX * 1e15; // Huge
    if (type == 3) return (double)rand() / RAND_MAX * DBL_MIN; // Subnormal boundary
    return (double)rand() / RAND_MAX * 100.0;
}

double rand_pos_double() {
    double x = rand_double();
    return x > 0 ? x : -x + 1e-9;
}

void test_ieee754_specials() {
    printf("--- IEEE 754 Special Values Gauntlet ---\n");
    double specials[] = {0.0, -0.0, INFINITY, -INFINITY, NAN, DBL_MAX, DBL_MIN, DBL_TRUE_MIN};
    int n = sizeof(specials)/sizeof(specials[0]);

    for(int i=0; i<n; i++) {
        double x = specials[i];
        int cls = ml_fp_classify(x);
        if (x == 0.0) CHECK(cls == 0 || cls == 1, "Zero classify");
        if (isinf(x)) CHECK(ml_isinf(x), "Inf check");
        if (isnan(x)) CHECK(ml_isnan(x), "NaN check");

        // Math functions shouldn't crash
        sine(x); cosine(x); tangent(x);
        exponential(x); logarithm(x > 0 ? x : 1.0);
        ml_sqrt(x >= 0 ? x : 0.0);
        ml_fast_rsqrt(x > 0 ? x : 1.0);
    }
    CHECK(1, "Survived IEEE specials");
}

void test_trig_identities() {
    printf("--- Trig Identities (10,000 iterations) ---\n");
    for(int i=0; i<10000; i++) {
        double x = rand_double();
        double s = sine(x);
        double c = cosine(x);

        CHECK_NEAR(s*s + c*c, 1.0, 1e-10, "Pythagorean identity");

        double s2 = sine(2.0 * x);
        double c2 = cosine(2.0 * x);
        CHECK_NEAR(s2, 2.0 * s * c, 1e-9, "sin(2x) identity");
        CHECK_NEAR(c2, c*c - s*s, 1e-9, "cos(2x) identity");

        CHECK_NEAR(sine(x + math_pi / 2.0), c, 1e-9, "sin(x+pi/2) == cos(x)");
    }
}

void test_exp_log_properties() {
    printf("--- Exp/Log Properties (5,000 iterations) ---\n");
    for(int i=0; i<5000; i++) {
        double x = rand_pos_double();
        double y = rand_pos_double();

        CHECK_NEAR(exponential(logarithm(x)), x, x * 1e-12, "exp(log(x)) == x");
        if (x < 700.0) {
            CHECK_NEAR(logarithm(exponential(x)), x, 1e-9, "log(exp(x)) == x");
        }

        if (x * y < 1e100) {
            CHECK_NEAR(logarithm(x * y), logarithm(x) + logarithm(y), 1e-9, "log(xy) == log(x)+log(y)");
        }

        double p = power(x, y);
        if (p < 1e100 && p > 1e-100) {
            CHECK_NEAR(logarithm(p), y * logarithm(x), 1e-8, "log(x^y) == y*log(x)");
        }
    }
}

void test_catastrophic_cancellation() {
    printf("--- Catastrophic Cancellation & Stability ---\n");

    double r1 = formula_pos(1.0, 1e8, 1.0);
    double r2 = formula_neg(1.0, 1e8, 1.0);
    CHECK_NEAR(r1, -1e-8, 1e-15, "Citardauq small root");
    CHECK_NEAR(r2, -1e8, 1.0, "Citardauq large root");

    double err;
    double s = ml_two_sum(1e16, 1.0, &err);
    CHECK_NEAR(s, 1e16, 1.0, "Two-Sum large");
    CHECK_NEAR(err, 1.0, 1e-15, "Two-Sum captured lost bit");

    double fma_res = ml_fma(1e16, 1.0, 1.0);
    CHECK_NEAR(fma_res, 1e16, 1.0, "FMA large");
}

void test_matrix_invariants() {
    printf("--- Matrix Invariants (500 iterations) ---\n");
    for(int i=0; i<500; i++) {
        mat3x3 A, B;
        for(int j=0; j<9; j++) {
            A.m[j] = (double)rand() / RAND_MAX * 10.0 - 5.0;
            B.m[j] = (double)rand() / RAND_MAX * 10.0 - 5.0;
        }

        double detA = mat3x3_det(A);
        double detB = mat3x3_det(B);
        mat3x3 AB = mat3x3_mul(A, B);
        double detAB = mat3x3_det(AB);

        if (fabs(detA) > 1e-5 && fabs(detB) > 1e-5) {
            CHECK_NEAR(detAB, detA * detB, 1e-4, "det(AB) == det(A)det(B)");
        }

        if (fabs(detA) > 1e-3) {
            mat3x3 invA = mat3x3_inverse(A);
            mat3x3 I = mat3x3_mul(A, invA);
            CHECK_NEAR(I.m[0], 1.0, 1e-6, "A*inv(A) == I [0]");
            CHECK_NEAR(I.m[4], 1.0, 1e-6, "A*inv(A) == I [4]");
            CHECK_NEAR(I.m[8], 1.0, 1e-6, "A*inv(A) == I [8]");
            CHECK_NEAR(I.m[1], 0.0, 1e-6, "A*inv(A) == I [1]");
        }
    }
}

void test_quaternion_algebra() {
    printf("--- Quaternion Algebra (1,000 iterations) ---\n");
    for(int i=0; i<1000; i++) {
        ml_quat q1 = {rand_double(), rand_double(), rand_double(), rand_double()};
        ml_quat q2 = {rand_double(), rand_double(), rand_double(), rand_double()};

        double n1 = q1.w*q1.w + q1.x*q1.x + q1.y*q1.y + q1.z*q1.z;
        double n2 = q2.w*q2.w + q2.x*q2.x + q2.y*q2.y + q2.z*q2.z;

        ml_quat q1q2 = ml_quat_mul(q1, q2);
        double n12 = q1q2.w*q1q2.w + q1q2.x*q1q2.x + q1q2.y*q1q2.y + q1q2.z*q1q2.z;

        CHECK_NEAR(n12, n1 * n2, 1e-6, "norm(q1*q2) == norm(q1)*norm(q2)");

        ml_quat conj1 = {q1.w, -q1.x, -q1.y, -q1.z};
        ml_quat q1_conj1 = ml_quat_mul(q1, conj1);
        CHECK_NEAR(q1_conj1.w, n1, 1e-6, "q*conj(q) == norm^2 (w)");
        CHECK_NEAR(q1_conj1.x, 0.0, 1e-6, "q*conj(q) == norm^2 (x)");
    }
}

void test_complex_identities() {
    printf("--- Complex Identities (1,000 iterations) ---\n");
    for(int i=0; i<1000; i++) {
        cplx z1 = {rand_double(), rand_double()};
        cplx z2 = {rand_double(), rand_double()};

        cplx prod = cplx_mul(z1, z2);
        double abs_prod = cplx_abs(prod);
        double abs1_abs2 = cplx_abs(z1) * cplx_abs(z2);

        CHECK_NEAR(abs_prod, abs1_abs2, 1e-7, "abs(z1*z2) == abs(z1)*abs(z2)");

        cplx euler = cplx_exponential((cplx){0.0, math_pi});
        CHECK_NEAR(euler.real, -1.0, 1e-9, "e^(i*pi) == -1 (real)");
        CHECK_NEAR(euler.imag, 0.0, 1e-9, "e^(i*pi) == 0 (imag)");
    }
}

void test_v10_tensor_solver() {
    printf("--- v10 Zero-Alloc Tensor Solver (100 iterations) ---\n");
    char scratchpad[1024 * 1024];
    ml_workspace_t ws = { scratchpad, sizeof(scratchpad), 0 };

    for(int iter=0; iter<100; iter++) {
        ml_workspace_reset(&ws);
        double A_data[16];
        double b_data[4] = {1.0, 2.0, 3.0, 4.0};
        double x_data[4] = {0};

        for(int i=0; i<4; i++) {
            double sum = 0;
            for(int j=0; j<4; j++) {
                if (i != j) {
                    A_data[i*4+j] = (double)rand() / RAND_MAX * 2.0 - 1.0;
                    sum += fabs(A_data[i*4+j]);
                }
            }
            A_data[i*4+i] = sum + 5.0;
        }

        ml_tensor_view_t A_view = ml_tensor_view(A_data, 4, 4);
        int status = ml_solve_v10(A_view, b_data, x_data, &ws);
        CHECK(status == 0, "Tensor solve status");

        for(int i=0; i<4; i++) {
            double sum = 0;
            for(int j=0; j<4; j++) sum += A_data[i*4+j] * x_data[j];
            CHECK_NEAR(sum, b_data[i], 1e-6, "Ax == b verification");
        }
    }
}

void test_simd_and_bare_metal() {
    printf("--- SIMD & Bare Metal Gauntlet ---\n");
    double in[1024] __attribute__((aligned(32)));
    double out_scalar[1024];
    double out_simd[1024] __attribute__((aligned(32)));

    for(int i=0; i<1024; i++) {
        in[i] = (double)rand() / RAND_MAX * 1000.0 + 1.0;
        out_scalar[i] = 1.0 / sqrt(in[i]);
    }

    for(int i=0; i<1024; i+=4) {
        ml_simd_batch_rsqrt(&in[i], &out_simd[i]);
    }
    for(int i=0; i<1024; i++) {
        CHECK_NEAR(out_simd[i], out_scalar[i], 1e-6, "SIMD batch rsqrt");
    }

    double out_avx[4] __attribute__((aligned(32)));
    __m256d v_in = _mm256_load_pd(in);
    __m256d v_out = ml_avx2_fast_rsqrt(v_in);
    _mm256_store_pd(out_avx, v_out);
    for(int i=0; i<4; i++) {
        CHECK_NEAR(out_avx[i], out_scalar[i], 1e-3, "AVX2 bare metal rsqrt");
    }
}

void test_fixed_point_cordic() {
    printf("--- Fixed-Point CORDIC vs Double CORDIC (1,000 iterations) ---\n");
    for(int i=0; i<1000; i++) {
        double angle = rand_double();
        while(angle > math_pi) angle -= 2.0*math_pi;
        while(angle < -math_pi) angle += 2.0*math_pi;

        ml_q16_16_t fixed_angle = (ml_q16_16_t)(angle * 65536.0);
        ml_q16_16_t f_sin, f_cos;
        ml_cordic_sincos_fixed(fixed_angle, &f_sin, &f_cos);

        double d_sin, d_cos;
        ml_cordic_sincos(angle, &d_sin, &d_cos);

        CHECK_NEAR((double)f_sin / 65536.0, d_sin, 1e-3, "Fixed vs Double CORDIC sin");
        CHECK_NEAR((double)f_cos / 65536.0, d_cos, 1e-3, "Fixed vs Double CORDIC cos");
    }
}

int main() {
    srand(time(NULL));
    printf("=========================================================\n");
    printf("   MATHLIB v10: GOD-MODE FUZZING GAUNTLET\n");
    printf("=========================================================\n\n");

    test_ieee754_specials();
    test_trig_identities();
    test_exp_log_properties();
    test_catastrophic_cancellation();
    test_matrix_invariants();
    test_quaternion_algebra();
    test_complex_identities();
    test_v10_tensor_solver();
    test_simd_and_bare_metal();
    test_fixed_point_cordic();

    printf("\n=========================================================\n");
    printf("GOD-MODE SUMMARY: %d passed, %d failed\n", passed, failed);
    printf("=========================================================\n");

    return failed > 0 ? 1 : 0;
}
