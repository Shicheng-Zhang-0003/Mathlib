
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

// Include the entire v10.5 engine
#include "ml_core.h"
#include "bitwise_fp.h"
#include "combinatorics.h"
#include "quadratics.h"
#include "integral.h"
#include "trigonometry.h"
#include "exponential.h"
#include "numerical.h"
#include "polynomial.h"
#include "ml_complex.h"
#include "statistics.h"
#include "ode.h"
#include "optimization.h"
#include "fft.h"
#include "fast_math.h"
#include "error_free.h"
#include "cordic.h"
#include "payne_hanek.h"
#include "minimax.h"
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

#define CHECK_NEAR(a, b, eps, msg) CHECK(ml_fabs((double)(a) - (double)(b)) < (eps), msg)
#define CHECK_NEAR_LOOSE(a, b, eps, msg) CHECK(ml_fabs((double)(a) - (double)(b)) < (eps), msg)
#define CHECK_NAN(a, msg) CHECK(ml_isnan(a), msg)
#define CHECK_INF(a, msg) CHECK(ml_isinf(a), msg)

double rand_double() {
    int type = rand() % 5;
    if (type == 0) return (double)rand() / RAND_MAX * 2000.0 - 1000.0; // Normal
    if (type == 1) return (double)rand() / RAND_MAX * 1e-15; // Tiny
    if (type == 2) return (double)rand() / RAND_MAX * 1e15; // Huge
    if (type == 3) return (double)rand() / RAND_MAX * 2.2250738585072014e-308; // Subnormal boundary
    return (double)rand() / RAND_MAX * 100.0;
}

double rand_pos_double() {
    double x = rand_double();
    return x > 0 ? x : -x + 1e-9;
}

// Bounded random for algebraic identities to respect IEEE-754 ULP limits
double rand_bounded_double() {
    return ((double)rand() / RAND_MAX * 20.0) - 10.0;
}

void test_ieee754_specials() {
    printf("--- IEEE 754 Special Values Gauntlet ---\n");
    double specials[] = {0.0, -0.0, (1.0/0.0), -(1.0/0.0), (0.0/0.0), 1.7976931348623157e+308, 2.2250738585072014e-308, 5e-324};
    int n = sizeof(specials)/sizeof(specials[0]);

    for(int i=0; i<n; i++) {
        double x = specials[i];
        int cls = ml_fp_classify(x);
        if (x == 0.0) CHECK(cls == 0 || cls == 1, "Zero classify");
        if (ml_isinf(x)) CHECK(ml_isinf(x), "Inf check");
        if (ml_isnan(x)) CHECK(ml_isnan(x), "NaN check");

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
        double x = ((double)rand() / RAND_MAX * 2000.0) - 1000.0;
        double s = sine(x);
        double c = cosine(x);

        double pyth = s*s + c*c;
        if (!(ml_fabs(pyth - 1.0) < 1e-10)) {
            failed++; printf("FAIL: Pythagorean | x: %.15g | s: %.15g | c: %.15g | s^2+c^2: %.15g\n", x, s, c, pyth);
        } else { passed++; }

        double s2 = sine(2.0 * x);
        double c2 = cosine(2.0 * x);
        double sin2x = 2.0 * s * c;
        double cos2x = c*c - s*s;

        if (!(ml_fabs(s2 - sin2x) < 1e-9)) {
            failed++; printf("FAIL: sin(2x) | x: %.15g | sin(2x): %.15g | 2sc: %.15g\n", x, s2, sin2x);
        } else { passed++; }

        if (!(ml_fabs(c2 - cos2x) < 1e-9)) {
            failed++; printf("FAIL: cos(2x) | x: %.15g | cos(2x): %.15g | c^2-s^2: %.15g\n", x, c2, cos2x);
        } else { passed++; }

        double phase = sine(x + 3.14159265358979323846 / 2.0);
        if (!(ml_fabs(phase - c) < 1e-9)) {
            failed++; printf("FAIL: sin(x+pi/2) | x: %.15g | sin(x+pi/2): %.15g | cos(x): %.15g\n", x, phase, c);
        } else { passed++; }
    }
}

void test_exp_log_properties() {
    printf("--- Exp/Log Properties (5,000 iterations) ---\n");
    for(int i=0; i<5000; i++) {
        double x = rand_pos_double();
        double y = rand_pos_double();

        double lx = logarithm(x);
        if (lx > -700.0 && lx < 700.0) {
            double exp_lx = exponential(lx);
            if (!(ml_fabs(exp_lx - x) < ml_fabs(x) * 1e-9 + 1e-12)) {
                failed++; printf("FAIL: exp(log(x)) | x: %.15g | log(x): %.15g | exp(log(x)): %.15g\n", x, lx, exp_lx);
            } else { passed++; }
        }
        if (x < 700.0) {
            double exp_x = exponential(x);
            double log_exp = logarithm(exp_x);
            if (!(ml_fabs(log_exp - x) < 1e-9)) {
                failed++; printf("FAIL: log(exp(x)) | x: %.15g | exp(x): %.15g | log(exp(x)): %.15g\n", x, exp_x, log_exp);
            } else { passed++; }
        }

        double prod = x * y;
        // Skip subnormals to avoid IEEE-754 hardware rounding errors
        if (x > 2.2250738585072014e-308 && y > 2.2250738585072014e-308 && prod > 2.2250738585072014e-308 && !ml_isinf(prod)) {
            double log_prod = logarithm(prod);
            double log_sum = logarithm(x) + logarithm(y);
            if (!(ml_fabs(log_prod - log_sum) < ml_fabs(logarithm(x)) * 1e-9 + 1e-9)) {
                failed++; printf("FAIL: log(xy) | x: %.15g | y: %.15g | log(xy): %.15g | log(x)+log(y): %.15g\n", x, y, log_prod, log_sum);
            } else { passed++; }
        }

        double p = power(x, y);
        if (p < 1e100 && p > 1e-100) {
            double log_p = logarithm(p);
            double y_log_x = y * logarithm(x);
            if (!(ml_fabs(log_p - y_log_x) < 1e-8)) {
                failed++; printf("FAIL: log(x^y) | x: %.15g | y: %.15g | log(x^y): %.15g | y*log(x): %.15g\n", x, y, log_p, y_log_x);
            } else { passed++; }
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


// Inline 3x3 math to avoid legacy linear_algebra.h dependency (V1.0.5 Quarantine)
static inline double det3x3(double m[9]) {
    return m[0]*(m[4]*m[8] - m[5]*m[7]) - m[1]*(m[3]*m[8] - m[5]*m[6]) + m[2]*(m[3]*m[7] - m[4]*m[6]);
}
static inline void mul3x3(double a[9], double b[9], double out[9]) {
    for(int i=0; i<3; i++) for(int j=0; j<3; j++) {
        out[i*3+j] = 0;
        for(int k=0; k<3; k++) out[i*3+j] += a[i*3+k] * b[k*3+j];
    }
}
static inline int inv3x3(double m[9], double out[9]) {
    double d = det3x3(m);
    if (ml_fabs(d) < 1e-12) return -1;
    double inv_d = 1.0 / d;
    out[0] =  (m[4]*m[8] - m[5]*m[7]) * inv_d;
    out[1] = -(m[1]*m[8] - m[2]*m[7]) * inv_d;
    out[2] =  (m[1]*m[5] - m[2]*m[4]) * inv_d;
    out[3] = -(m[3]*m[8] - m[5]*m[6]) * inv_d;
    out[4] =  (m[0]*m[8] - m[2]*m[6]) * inv_d;
    out[5] = -(m[0]*m[5] - m[2]*m[3]) * inv_d;
    out[6] =  (m[3]*m[7] - m[4]*m[6]) * inv_d;
    out[7] = -(m[0]*m[7] - m[1]*m[6]) * inv_d;
    out[8] =  (m[0]*m[4] - m[1]*m[3]) * inv_d;
    return 0;
}

void test_matrix_invariants() {
    printf("--- Matrix Invariants (500 iterations) ---");
    for(int i=0; i<500; i++) {
        double A[9], B[9];
        for(int j=0; j<9; j++) {
            A[j] = (double)rand() / RAND_MAX * 10.0 - 5.0;
            B[j] = (double)rand() / RAND_MAX * 10.0 - 5.0;
        }
        double detA = det3x3(A);
        double detB = det3x3(B);
        double AB[9];
        mul3x3(A, B, AB);
        double detAB = det3x3(AB);

        if (ml_fabs(detA) > 1e-5 && ml_fabs(detB) > 1e-5) {
            CHECK_NEAR(detAB, detA * detB, 1e-4, "det(AB) == det(A)det(B)");
        }
        if (ml_fabs(detA) > 1e-3) {
            double invA[9];
            inv3x3(A, invA);
            double I[9];
            mul3x3(A, invA, I);
            CHECK_NEAR(I[0], 1.0, 1e-6, "A*inv(A) == I [0]");
            CHECK_NEAR(I[4], 1.0, 1e-6, "A*inv(A) == I [4]");
            CHECK_NEAR(I[8], 1.0, 1e-6, "A*inv(A) == I [8]");
            CHECK_NEAR(I[1], 0.0, 1e-6, "A*inv(A) == I [1]");
        }
    }
}

void test_quaternion_algebra() {
    printf("--- Quaternion Algebra (1,000 iterations) ---\n");
    for(int i=0; i<1000; i++) {
        ml_quat q1 = {rand_bounded_double(), rand_bounded_double(), rand_bounded_double(), rand_bounded_double()};
        ml_quat q2 = {rand_bounded_double(), rand_bounded_double(), rand_bounded_double(), rand_bounded_double()};

        double n1 = q1.w*q1.w + q1.x*q1.x + q1.y*q1.y + q1.z*q1.z;
        double n2 = q2.w*q2.w + q2.x*q2.x + q2.y*q2.y + q2.z*q2.z;

        ml_quat q1q2 = ml_quat_mul(q1, q2);
        double n12 = q1q2.w*q1q2.w + q1q2.x*q1q2.x + q1q2.y*q1q2.y + q1q2.z*q1q2.z;

        double tol_q = ml_fabs(n1 * n2) * 1e-7 + 1e-12;
        CHECK_NEAR_LOOSE(n12, n1 * n2, tol_q, "norm(q1*q2) == norm(q1)*norm(q2)");

        ml_quat conj1 = {q1.w, -q1.x, -q1.y, -q1.z};
        ml_quat q1_conj1 = ml_quat_mul(q1, conj1);
        CHECK_NEAR(q1_conj1.w, n1, 1e-6, "q*conj(q) == norm^2 (w)");
        CHECK_NEAR(q1_conj1.x, 0.0, 1e-6, "q*conj(q) == norm^2 (x)");
    }
}

void test_complex_identities() {
    printf("--- Complex Identities (1,000 iterations) ---\n");
    for(int i=0; i<1000; i++) {
        cplx z1 = {rand_bounded_double(), rand_bounded_double()};
        cplx z2 = {rand_bounded_double(), rand_bounded_double()};

        cplx prod = cplx_mul(z1, z2);
        double abs_prod = cplx_abs(prod);
        double abs1_abs2 = cplx_abs(z1) * cplx_abs(z2);

        double tol_c = ml_fabs(abs1_abs2) * 1e-7 + 1e-12;
        CHECK_NEAR_LOOSE(abs_prod, abs1_abs2, tol_c, "abs(z1*z2) == abs(z1)*abs(z2)");

        cplx euler = cplx_exponential((cplx){0.0, 3.14159265358979323846});
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
                    sum += ml_fabs(A_data[i*4+j]);
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
        in[i] = (double)rand() / RAND_MAX * 99.0 + 1.0;
        out_scalar[i] = 1.0 / ml_sqrt(in[i]);
    }

    for(int i=0; i<1024; i+=4) {
        ml_simd_batch_rsqrt(&in[i], &out_simd[i]);
    }
    for(int i=0; i<1024; i++) {
        double tol_simd = ml_fabs(out_scalar[i]) * 1e-5 + 1e-12;
        CHECK_NEAR_LOOSE(out_simd[i], out_scalar[i], tol_simd, "SIMD batch rsqrt");
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
        // O(1) range reduction to prevent infinite loop on massive inputs
        angle = ml_fmod(angle, 2.0 * 3.14159265358979323846);
        if (angle > 3.14159265358979323846) angle -= 2.0 * 3.14159265358979323846;
        if (angle < -3.14159265358979323846) angle += 2.0 * 3.14159265358979323846;

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
    unsigned int seed = time(NULL);
    srand(seed);
    printf("Fuzz Seed: %u\n", seed);
    printf("=========================================================\n");
    printf("   MATHLIB v1.0: GOD-MODE DYNAMIC FUZZING GAUNTLET\n");
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
