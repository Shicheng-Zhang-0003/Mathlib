
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

// Fuzz drivers below exercise the staged ml_* kernels directly.
#include "ml_core.h"
#include "bitwise_fp.h"
#include "ml_combinatorics.h"
#include "ml_quadratics.h"
#include "ml_integral.h"
#include "ml_trig.h"
#include "ml_exp_log.h"
#include "ml_numerical.h"
#include "ml_polynomial.h"
#include "ml_complex.h"
#include "ml_statistics.h"
#include "ml_ode.h"
#include "ml_optimization.h"
#include "fft.h"
#include "fast_math.h"
#include "internal/error_free.h"
#include "internal/cordic.h"
#include "internal/payne_hanek.h"
#include "internal/minimax.h"
#include "ieee754.h"
#include "simd.h"
#include "ml_quaternion.h"
#include "ml_tensor.h"
#include "ml_linalg.h"
#include "ml_types.h"
#include "ml_fixed_point.h"
#include "ml_numbertheory.h"
#include "ml_contest.h"
#include "ml_geometry.h"
#include "ml_graph.h"
#include "ml_modeling.h"
#include "ml_measure.h"
#include "ml_funcanal.h"
#include "ml_topology.h"
#include "ml_diffgeom.h"
#include "ml_algebraic.h"
#include "ml_stats_inv.h"
#include "ml_spectral.h"
#include "ml_ode.h"
#include "ml_pde.h"
#include "ml_harmonic.h"
#include "ml_calculus.h"
#include "ml_transforms.h"
#include "ml_info.h"
#include "ml_mcmc.h"
#include "ml_manifold.h"
#include "ml_quadratics.h"
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
    int type = rand() % 10;
    if (type == 0) return (double)rand() / RAND_MAX * 2000.0 - 1000.0; // Normal [-1000, 1000]
    if (type == 1) return (double)rand() / RAND_MAX * 2e-15 - 1e-15;  // Tiny [-1e-15, 1e-15]
    if (type == 2) return (double)rand() / RAND_MAX * 2e15 - 1e15;    // Huge [-1e15, 1e15]
    if (type == 3) return (double)rand() / RAND_MAX * 4.45e-308 - 2.225e-308; // Subnormal boundary
    if (type == 4) return 0.0;                                        // Exact zero
    if (type == 5) return -0.0;                                       // Negative zero
    if (type == 6) return ml_make_inf(0);                                  // Infinity
    if (type == 7) return -ml_make_inf(0);                                 // Negative Infinity
    if (type == 8) return ml_make_nan();                                  // NaN
    return (double)rand() / RAND_MAX * 200.0 - 100.0;                 // Standard [-100, 100]
}

double rand_pos_double() {
    double x = rand_double();
    return x > 0 ? x : -x + 1e-9;
}

// Bounded random for algebraic identities to respect IEEE-754 ULP limits
double rand_bounded_double() {
    return ((double)rand() / RAND_MAX * 20.0) - 10.0;
}

/* MATHLIB_CLOSURE_P2_LOG_EXP_ROUNDTRIP_HELPER */
static uint64_t ml_fuzz_ulp_distance(double a, double b) {
    uint64_t ia, ib;

    memcpy(&ia, &a, sizeof(uint64_t));
    memcpy(&ib, &b, sizeof(uint64_t));

    if (ia >> 63) ia = 0x8000000000000000ULL - ia;
    if (ib >> 63) ib = 0x8000000000000000ULL - ib;

    return ia > ib ? ia - ib : ib - ia;
}

void test_ieee754_specials() {
    printf("--- IEEE 754 Special Values Gauntlet ---\n");
    double specials[] = {0.0, -0.0, (ml_make_inf(0)), -(ml_make_inf(0)), (ml_make_nan()), 1.7976931348623157e+308, 2.2250738585072014e-308, 5e-324};
    int n = sizeof(specials)/sizeof(specials[0]);

    for(int i=0; i<n; i++) {
        double x = specials[i];
        int cls = ml_fp_classify(x);
        if (x == 0.0) CHECK(cls == 0 || cls == 1, "Zero classify");
        if (ml_isinf(x)) CHECK(ml_isinf(x), "Inf check");
        if (ml_isnan(x)) CHECK(ml_isnan(x), "NaN check");

        ml_sin(x); ml_cos(x); ml_tan(x);
        ml_exp(x); ml_log(x > 0 ? x : 1.0);
        ml_sqrt(x >= 0 ? x : 0.0);
        ml_fast_rsqrt(x > 0 ? x : 1.0);
    }
    CHECK(1, "Survived IEEE specials");
}

void test_trig_identities() {
    printf("--- Trig Identities (10,000 iterations) ---\n");
    for(int i=0; i<10000; i++) {
        double x = ((double)rand() / RAND_MAX * 2000.0) - 1000.0;
        double s = ml_sin(x);
        double c = ml_cos(x);

        double pyth = s*s + c*c;
        if (!(ml_fabs(pyth - 1.0) < 1e-14)) {
            failed++; printf("FAIL: Pythagorean | x: %.15g | s: %.15g | c: %.15g | s^2+c^2: %.15g\n", x, s, c, pyth);
        } else { passed++; }

        double s2 = ml_sin(2.0 * x);
        double c2 = ml_cos(2.0 * x);
        double sin2x = 2.0 * s * c;
        double cos2x = c*c - s*s;

        if (!(ml_fabs(s2 - sin2x) < 1e-13)) {
            failed++; printf("FAIL: sin(2x) | x: %.15g | sin(2x): %.15g | 2sc: %.15g\n", x, s2, sin2x);
        } else { passed++; }

        if (!(ml_fabs(c2 - cos2x) < 1e-13)) {
            failed++; printf("FAIL: cos(2x) | x: %.15g | cos(2x): %.15g | c^2-s^2: %.15g\n", x, c2, cos2x);
        } else { passed++; }

        double phase = ml_sin(x + ML_PI / 2.0);
        if (!(ml_fabs(phase - c) < 1e-13)) {
            failed++; printf("FAIL: sin(x+pi/2) | x: %.15g | sin(x+pi/2): %.15g | cos(x): %.15g\n", x, phase, c);
        } else { passed++; }
    }
}

void test_exp_log_properties() {
    printf("--- Exp/Log Properties (5,000 iterations) ---\n");
    for(int i=0; i<5000; i++) {
        double x = rand_pos_double();
        double y = rand_pos_double();

        double lx = ml_log(x);
        if (lx > -700.0 && lx < 700.0) {
            double exp_lx = ml_exp(lx);
            if (!(ml_fabs(exp_lx - x) < ml_fabs(x) * 1e-14 + 1e-15)) {
                failed++; printf("FAIL: exp(log(x)) | x: %.15g | log(x): %.15g | exp(log(x)): %.15g\n", x, lx, exp_lx);
            } else { passed++; }
        }
        if (x < 700.0) {
    double exp_x = ml_exp(x);
    double log_exp = ml_log(exp_x);
    double diff = ml_fabs(log_exp - x);
    uint64_t ulps = ml_fuzz_ulp_distance(log_exp, x);

    /*
     * MATHLIB_CLOSURE_P2_LOG_EXP_ROUNDTRIP_TOLERANCE
     *
     * The old absolute tolerance 1e-13 is smaller than one double ULP
     * for x in [512, 1024), where one ULP is approximately 1.14e-13.
     *
     * Therefore a roundtrip error of one ULP could fail the old test
     * even though the result is numerically sane.
     *
     * Accept either:
     *   - very small absolute error, or
     *   - <= 8 ULP roundtrip error.
     */
    if (ml_isnan(log_exp) || ml_isinf(log_exp) || !(diff <= 1e-13 || ulps <= 8)) {
        failed++;
        printf("FAIL: log(exp(x)) | x: %.17g | exp(x): %.17g | log(exp(x)): %.17g | diff: %.17g | ulps: %llu\n",
               x, exp_x, log_exp, diff, (unsigned long long)ulps);
    } else {
        passed++;
    }
}

double prod = x * y;
        // Skip subnormals to avoid IEEE-754 hardware rounding errors
        if (x > 2.2250738585072014e-308 && y > 2.2250738585072014e-308 && prod > 2.2250738585072014e-308 && !ml_isinf(prod)) {
            double log_prod = ml_log(prod);
            double log_sum = ml_log(x) + ml_log(y);
            if (!(ml_fabs(log_prod - log_sum) < ml_fabs(ml_log(x)) * 1e-13 + 1e-14)) {
                failed++; printf("FAIL: log(xy) | x: %.15g | y: %.15g | log(xy): %.15g | log(x)+log(y): %.15g\n", x, y, log_prod, log_sum);
            } else { passed++; }
        }

        double p = ml_pow(x, y);
        if (p < 1e100 && p > 1e-100) {
            double log_p = ml_log(p);
            double y_log_x = y * ml_log(x);
            if (!(ml_fabs(log_p - y_log_x) < 1e-13)) {
                failed++; printf("FAIL: log(x^y) | x: %.15g | y: %.15g | log(x^y): %.15g | y*log(x): %.15g\n", x, y, log_p, y_log_x);
            } else { passed++; }
        }
    }
}
void test_catastrophic_cancellation() {
    printf("--- Catastrophic Cancellation & Stability ---\n");

    double r1 = ml_formula_pos(1.0, 1e8, 1.0);
    double r2 = ml_formula_neg(1.0, 1e8, 1.0);
    CHECK_NEAR(r1, -1e-8, 1e-15, "Citardauq small root");
    CHECK_NEAR(r2, -1e8, 1.0, "Citardauq large root");

    double err;
    double s = ml_two_sum(1e16, 1.0, &err);
    CHECK_NEAR(s, 1e16, 1.0, "Two-Sum large");
    CHECK_NEAR(err, 1.0, 1e-15, "Two-Sum captured lost bit");

    double fma_res = ml_fma_soft(1e16, 1.0, 1.0);
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

        cplx prod = ml_cplx_mul(z1, z2);
        double abs_prod = ml_cplx_abs(prod);
        double abs1_abs2 = ml_cplx_abs(z1) * ml_cplx_abs(z2);

        double tol_c = ml_fabs(abs1_abs2) * 1e-7 + 1e-12;
        CHECK_NEAR_LOOSE(abs_prod, abs1_abs2, tol_c, "abs(z1*z2) == abs(z1)*abs(z2)");

        cplx euler = ml_cplx_exponential((cplx){0.0, ML_PI});
        CHECK_NEAR(euler.real, -1.0, 1e-9, "e^(i*pi) == -1 (real)");
        CHECK_NEAR(euler.imag, 0.0, 1e-9, "e^(i*pi) == 0 (imag)");
    }
}

void test_v10_tensor_solver() {
    printf("--- v10 Zero-Alloc Tensor Solver (100 iterations) ---\n");
    char scratchpad[1024 * 1024];
    ml_workspace_t ws;
    ml_workspace_init(&ws, scratchpad, sizeof(scratchpad));

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
        int status = ml_solve(A_view, b_data, x_data, &ws);
        CHECK(status == ML_SUCCESS, "Tensor solve status");

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

#if defined(__AVX__) || defined(__AVX2__)
    double out_avx[4] __attribute__((aligned(32)));
    __m256d v_in = _mm256_load_pd(in);
    __m256d v_out = ml_avx2_fast_rsqrt(v_in);
    _mm256_store_pd(out_avx, v_out);
    for(int i=0; i<4; i++) {
        CHECK_NEAR(out_avx[i], out_scalar[i], 1e-3, "AVX2 bare metal rsqrt");
    }
#else
    printf("  [SKIP] AVX2 bare metal rsqrt (Non-AVX build)\n");
#endif
}

void test_fixed_point_cordic() {
    /* MATHLIB_CLOSURE_P2_P0_5_FUZZ_CORDIC_NONFINITE_GUARD */
    printf("--- Fixed-Point CORDIC vs Double CORDIC (1,000 iterations) ---\n");

    for(int i=0; i<1000; i++) {
        double angle = rand_double();

        /*
         * Do NOT cast NaN or Infinity to integer.
         *
         * In C, converting a non-finite floating-point value to an integer
         * type is undefined behavior.
         *
         * Instead, verify that the double CORDIC path correctly returns NaN.
         */
        if (!ml_isfinite(angle)) {
            double d_sin, d_cos;
            ml_cordic_sincos(angle, &d_sin, &d_cos);

            CHECK(ml_isnan(d_sin) && ml_isnan(d_cos),
                  "Double CORDIC non-finite input produces NaN");

            continue;
        }

        /* O(1) range reduction to prevent infinite loop on massive inputs */
        angle = ml_fmod(angle, 2.0 * ML_PI);

        if (!ml_isfinite(angle)) {
            CHECK(0, "ml_fmod(finite, 2*pi) produced non-finite result");
            continue;
        }

        if (angle > ML_PI) angle -= 2.0 * ML_PI;
        if (angle < -ML_PI) angle += 2.0 * ML_PI;

        ml_q16_16_t fixed_angle = (ml_q16_16_t)(angle * 65536.0);

        ml_q16_16_t f_sin, f_cos;
        ml_cordic_sincos_fixed(fixed_angle, &f_sin, &f_cos);

        double d_sin, d_cos;
        ml_cordic_sincos(angle, &d_sin, &d_cos);

        CHECK_NEAR((double)f_sin / 65536.0, d_sin, 1e-3, "Fixed vs Double CORDIC sin");
        CHECK_NEAR((double)f_cos / 65536.0, d_cos, 1e-3, "Fixed vs Double CORDIC cos");
    }
}


/* ============================================================================
 * GOD-MODE PART 2 (fuzz upgrade): invariant + differential fuzz for every
 * module the original gauntlet never touched (Batch-1, contest, geometry,
 * graph, modeling, doctorate). Discipline differs from Part 1 on purpose:
 * - Self-contained xorshift64* stream (fx_*): counts are host-stable
 *   (Part 1 keeps rand(), untouched, so its stream and counts are unchanged).
 * - Fixed checks per iteration: inputs are constructed valid (bounded
 *   rejection sampling with a deterministic fallback), so the total is
 *   seed-independent; the gate stays "0 failed".
 * Each section header states its exact check budget.
 * ========================================================================== */

static uint64_t fx_state = 0x9E3779B97F4A7C15ULL;

static uint64_t fx_next(void) {
    uint64_t x = fx_state;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    fx_state = x;
    return x * 0x2545F4914F6CDD1DULL;
}

static double fx_unit(void) {
    return (double)(fx_next() >> 11) * (1.0 / 9007199254740992.0);
}

static uint64_t fx_below(uint64_t n) {
    return n == 0 ? 0 : fx_next() % n;
}

static double fx_range(double lo, double hi) {
    return lo + (hi - lo) * fx_unit();
}

static double fx_ode_exp(double t, double y) {
    (void)t;
    return y;
}

/* S1: contest identities, 800 iters x 7 = 5600. */
void test_fx_contest_identities(void) {
    static const int ones8[8] = {1, 1, 1, 1, 1, 1, 1, 1};
    int it;
    printf("--- Contest Identities (800 iterations) ---\n");
    for (it = 0; it < 800; it++) {
        int nn = 1 + (int)fx_below(8);
        int kk;
        uint64_t row = 0;
        for (kk = 1; kk <= nn; kk++) row += ml_narayana(nn, kk);
        CHECK(row == ml_catalan(nn), "fx narayana row == catalan");
        {
            int single[1];
            single[0] = nn;
            CHECK(ml_hook_dimension(single, 1) == 1, "fx hook single row");
            CHECK(ml_sn_class_size(ones8, nn) == 1, "fx class identity");
        }
        {
            /* Lucas mod 2 vs direct C(n,k) parity, n <= 30 exact. */
            int n2 = (int)fx_below(31);
            int k2 = (n2 == 0) ? 0 : (int)fx_below((uint64_t)n2 + 1U);
            uint64_t cnk = 1;
            int ii;
            for (ii = 1; ii <= k2; ii++) cnk = cnk * (uint64_t)(n2 - k2 + ii) / (uint64_t)ii;
            CHECK((ml_binom_odd((uint64_t)n2, (uint64_t)k2) != 0) == ((cnk & 1ULL) != 0), "fx lucas parity");
        }
        {
            /* Josephus vs naive simulation, m <= 20. */
            int m = 1 + (int)fx_below(20);
            int k = 1 + (int)fx_below((uint64_t)m);
            int alive[20];
            int rem;
            int idx = 0;
            int r;
            for (r = 0; r < m; r++) alive[r] = r;
            for (rem = m; rem > 1; rem--) {
                int q;
                idx = (idx + k - 1) % rem;
                for (q = idx; q < rem - 1; q++) alive[q] = alive[q + 1];
            }
            CHECK(ml_josephus((uint64_t)m, (uint64_t)k) == (uint64_t)alive[0], "fx josephus sim");
        }
        {
            /* Hanoi recurrence + Lucas/Fibonacci relation L(n) = F(n-1)+F(n+1). */
            int n3 = 1 + (int)fx_below(30);
            int n4 = 2 + (int)fx_below(19);
            uint64_t fa = 0, fb = 0, fc = 0, fd = 0;
            ml_fib_pair((long long)n4 - 1, &fa, &fb);
            ml_fib_pair((long long)n4 + 1, &fc, &fd);
            (void)fb;
            (void)fd;
            CHECK(ml_hanoi_moves(n3) == 2 * ml_hanoi_moves(n3 - 1) + 1, "fx hanoi recurrence");
            CHECK(ml_lucas_number(n4) == fa + fc, "fx lucas-fib relation");
        }
    }
}

/* S2: olympiad NT differentials, 600 iters x 7 = 4200. */
void test_fx_olympiad_nt(void) {
    static const uint64_t primes8[8] = {7, 11, 13, 17, 19, 23, 29, 31};
    static const uint64_t oddp[3] = {3, 5, 7};
    int it;
    printf("--- Olympiad NT Differentials (600 iterations) ---\n");
    for (it = 0; it < 600; it++) {
        {
            /* floor_sum vs naive, tiny exact. */
            uint64_t n = fx_below(30);
            uint64_t m = 1 + fx_below(50);
            uint64_t a = fx_below(50);
            uint64_t b = fx_below(50);
            uint64_t naive = 0;
            uint64_t i;
            for (i = 0; i < n; i++) naive += (a * i + b) / m;
            CHECK(ml_floor_sum_u64(n, m, a, b) == naive, "fx floor_sum naive");
        }
        if ((it & 1) == 0) {
            /* Odd-p LTE vs direct valuation, a,b <= 55, n <= 6 exact. */
            uint64_t p = oddp[fx_below(3)];
            uint64_t b = 1 + fx_below(20);
            uint64_t t = 1 + fx_below(5);
            uint64_t a = b + p * t;
            uint64_t n = 1 + fx_below(6);
            uint64_t pa = 1, pb = 1;
            uint64_t i;
            uint64_t diff;
            uint64_t got;
            uint64_t want = 0;
            if (a % p == 0) {
                a = 10;
                b = 4;
                p = 3;
                n = 6;
            }
            for (i = 0; i < n; i++) {
                pa *= a;
                pb *= b;
            }
            diff = (pa > pb) ? pa - pb : pb - pa;
            while (diff % p == 0) {
                diff /= p;
                want++;
            }
            got = ml_lte_vp(a, b, n, p);
            CHECK(got == want, "fx LTE direct");
        } else {
            /* 2-adic LTE, a,b odd, n even (a==b falls back: v2(0) undefined). */
            uint64_t a = 2 * fx_below(25) + 1;
            uint64_t b = 2 * fx_below(25) + 1;
            uint64_t n = 2 * (1 + fx_below(3));
            uint64_t pa = 1, pb = 1;
            uint64_t i;
            uint64_t diff;
            uint64_t want = 0;
            if (a == b) {
                a = 5;
                b = 3;
            }
            for (i = 0; i < n; i++) {
                pa *= a;
                pb *= b;
            }
            diff = (pa > pb) ? pa - pb : pb - pa;
            while ((diff & 1ULL) == 0 && diff != 0) {
                diff /= 2;
                want++;
            }
            CHECK(ml_lte_vp(a, b, n, 2) == want, "fx LTE 2-adic");
        }
        {
            /* Pell verify by substitution, D <= 30 non-square. */
            uint64_t D = 2 + fx_below(29);
            uint64_t s = 2;
            uint64_t x = 0, y = 0;
            while (s * s < D) s++;
            if (s * s == D) D = 2;
            if (ml_pell_fund(D, &x, &y) == ML_SUCCESS) {
                CHECK(x * x - D * y * y == 1, "fx pell verify");
            } else {
                CHECK(0, "fx pell solves D<=30");
            }
        }
        {
            /* Tonelli: least root or certified non-residue. */
            uint64_t p = primes8[fx_below(8)];
            uint64_t a = 1 + fx_below(p - 1);
            uint64_t r = ml_modsqrt_tonelli(a, p);
            if (r == UINT64_MAX) {
                CHECK(ml_modpow(a, (p - 1) / 2, p) == p - 1, "fx tonelli non-residue");
            } else {
                CHECK((r * r) % p == a && r <= p - r, "fx tonelli verify");
            }
        }
        {
            /* CF convergent quality |D-(p/q)^2| < 2/q^2. */
            uint64_t D = 2 + fx_below(49);
            uint64_t s = 2;
            uint64_t num = 0, den = 0;
            int k = (int)fx_below(7);
            while (s * s < D) s++;
            if (s * s == D) D = 2;
            CHECK(ml_cf_convergent_sqrt(D, k, &num, &den) == ML_SUCCESS, "fx cf status");
            {
                /* Convergents: |D q^2 - p^2| <= 2*floor(sqrt(D)) + 1, exact
                 * integer check (the 2/q^2-style bound fails at k = 0). */
                uint64_t s2 = 2;
                uint64_t lhs;
                uint64_t bound;
                while ((s2 + 1) * (s2 + 1) <= D) s2++;
                bound = 2 * s2 + 1;
                lhs = (D * den * den > num * num) ? D * den * den - num * num
                                                  : num * num - D * den * den;
                CHECK(lhs <= bound, "fx cf quality");
            }
        }
        {
            /* vp(n!) vs naive accumulation, n <= 60. */
            uint64_t n = fx_below(61);
            uint64_t p = oddp[fx_below(3)];
            uint64_t naive = 0;
            uint64_t j;
            for (j = 2; j <= n; j++) {
                uint64_t t2 = j;
                while (t2 % p == 0) {
                    t2 /= p;
                    naive++;
                }
            }
            CHECK(ml_vp_factorial((long long)n, p) == naive, "fx vp_fact naive");
        }
    }
}

/* S3: graph invariants, 500 iters x 9 = 4500. */
void test_fx_graph_invariants(void) {
    int it;
    printf("--- Graph Invariants (500 iterations) ---\n");
    for (it = 0; it < 500; it++) {
        int n = 2 + (int)fx_below(7);
        int eu[16], ev[16];
        int m = 0;
        int v;
        int e2;
        for (v = 1; v < n; v++) {
            eu[m] = (int)fx_below((uint64_t)v);
            ev[m] = v;
            m++;
        }
        for (e2 = 0; e2 < 6; e2++) {
            int a = (int)fx_below((uint64_t)n);
            int b = (int)fx_below((uint64_t)n);
            if (a == b) continue;
            eu[m] = a;
            ev[m] = b;
            m++;
        }
        CHECK(ml_components(n, eu, ev, m) == 1, "fx connected");
        {
            double dist[8];
            int ok = 1;
            int j;
            ml_bfs_dist(n, eu, ev, m, 0, dist);
            CHECK(dist[0] == 0.0, "fx bfs src");
            for (j = 0; j < n; j++) {
                if (!(dist[j] >= 0.0)) {
                    ok = 0;
                    break;
                }
            }
            CHECK(ok, "fx bfs reachable");
            {
                /* Dijkstra on unit weights == BFS hops, exactly. */
                double w[16];
                double dd[8];
                int same = 1;
                for (j = 0; j < m; j++) w[j] = 1.0;
                ml_dijkstra(n, eu, ev, w, m, 0, dd);
                for (j = 0; j < n; j++) {
                    if (dd[j] != dist[j]) {
                        same = 0;
                        break;
                    }
                }
                CHECK(same, "fx dijkstra==bfs unit");
            }
            {
                /* Kruskal on unit weights == n-1 exactly. */
                double w[16];
                double tot = 0.0;
                int j2;
                for (j2 = 0; j2 < m; j2++) w[j2] = 1.0;
                CHECK(ml_kruskal_mst(n, eu, ev, w, m, &tot) == ML_SUCCESS, "fx kruskal st");
                CHECK(tot == (double)(n - 1), "fx kruskal unit n-1");
            }
        }
        {
            /* Random DAG: permutation + forward edges. */
            int n2 = 2 + (int)fx_below(6);
            int perm[8];
            int du[28], dv[28];
            int m2 = 0;
            int a2;
            int b2;
            int out[8];
            int pos[8];
            int cnt;
            int ok = 1;
            for (a2 = 0; a2 < n2; a2++) perm[a2] = a2;
            for (a2 = n2 - 1; a2 > 0; a2--) {
                int bb = (int)fx_below((uint64_t)(a2 + 1));
                int tmp = perm[a2];
                perm[a2] = perm[bb];
                perm[bb] = tmp;
            }
            for (a2 = 0; a2 < n2; a2++) {
                for (b2 = a2 + 1; b2 < n2; b2++) {
                    if (fx_unit() < 0.4 && m2 < 28) {
                        du[m2] = perm[a2];
                        dv[m2] = perm[b2];
                        m2++;
                    }
                }
            }
            cnt = ml_topo_sort(n2, du, dv, m2, out);
            CHECK(cnt == n2, "fx topo count");
            for (a2 = 0; a2 < n2; a2++) pos[out[a2]] = a2;
            for (a2 = 0; a2 < m2; a2++) {
                if (!(pos[du[a2]] < pos[dv[a2]])) {
                    ok = 0;
                    break;
                }
            }
            CHECK(ok, "fx topo valid");
        }
        {
            /* Cycle parity: even <-> bipartite. */
            int L = 3 + (int)fx_below(6);
            int cu[8], cv[8];
            int j;
            for (j = 0; j < L; j++) {
                cu[j] = j;
                cv[j] = (j + 1) % L;
            }
            CHECK(ml_is_bipartite(L, cu, cv, L) == (L % 2 == 0 ? 1 : 0), "fx cycle parity");
        }
    }
}

/* S4: modeling closed forms, 500 iters x 9 = 4500. */
void test_fx_modeling_forms(void) {
    int it;
    printf("--- Modeling Closed Forms (500 iterations) ---\n");
    for (it = 0; it < 500; it++) {
        double P[4];
        double pi[2] = {0.5, 0.5};
        double o[2] = {0.0, 0.0};
        double r0;
        double s;
        double s2 = 0.0;
        int j;
        for (j = 0; j < 4; j++) P[j] = 0.05 + fx_unit();
        {
            double rs = P[0] + P[1];
            P[0] /= rs;
            P[1] /= rs;
            rs = P[2] + P[3];
            P[2] /= rs;
            P[3] /= rs;
        }
        CHECK(ml_markov_stationary(P, pi, 2, 10000, 1e-12) == ML_SUCCESS, "fx stationary st");
        s = pi[0] + pi[1];
        CHECK(ml_fabs(s - 1.0) < 1e-9, "fx stationary sums to 1");
        {
            double r1 = ml_fabs((pi[0] * P[0] + pi[1] * P[2]) - pi[0]);
            double r2 = ml_fabs((pi[0] * P[1] + pi[1] * P[3]) - pi[1]);
            CHECK(r1 < 1e-9 && r2 < 1e-9 && pi[0] >= 0.0 && pi[1] >= 0.0, "fx stationary residual");
        }
        CHECK(ml_markov_step(P, pi, o, 2) == ML_SUCCESS, "fx step st");
        CHECK(ml_fabs((o[0] + o[1]) - 1.0) < 1e-12, "fx step preserves sum");
        {
            /* Little's law L = lambda W + SIR final-size residual. */
            double lam = 0.5 + 9.5 * fx_unit();
            double mu = lam + 0.5 + 9.5 * fx_unit();
            double L = ml_mm1_L(lam, mu);
            double W = ml_mm1_W(lam, mu);
            CHECK(ml_fabs(L - lam * W) < 1e-9 * (1.0 + L), "fx little law");
            r0 = 1.1 + 8.9 * fx_unit();
            {
                double z = ml_sir_final_size(r0);
                double res = ml_fabs(z - (1.0 - ml_exp(-r0 * z)));
                CHECK(res < 1e-12, "fx sir residual");
            }
        }
        {
            /* Conjugate exactness: a1-a0 == succ. */
            uint64_t succ = fx_below(50);
            uint64_t tri = succ + fx_below(50);
            double a1 = 0, b1 = 0;
            CHECK(ml_beta_binomial_post(2, 2, succ, tri, &a1, &b1) == ML_SUCCESS, "fx conjugate st");
            CHECK(a1 - 2.0 == (double)succ, "fx conjugate exact");
        }
        (void)s2;
    }
}

/* S5: geometry predicates, 600 iters x 9 = 5400. */
void test_fx_geometry_predicates(void) {
    int it;
    printf("--- Geometry Predicates (600 iterations) ---\n");
    for (it = 0; it < 600; it++) {
        double x1 = fx_range(-10.0, 10.0);
        double y1 = fx_range(-10.0, 10.0);
        double x2 = fx_range(-10.0, 10.0);
        double y2 = fx_range(-10.0, 10.0);
        double x3 = fx_range(-10.0, 10.0);
        double y3 = fx_range(-10.0, 10.0);
        double area;
        double a;
        double b;
        double c;
        double R;
        double r;
        /* Fallback fixed triangle keeps every check unconditional. */
        double cross = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1);
        if (ml_fabs(cross) < 1e-6) {
            x1 = 0;
            y1 = 0;
            x2 = 4;
            y2 = 0;
            x3 = 0;
            y3 = 3;
        }
        area = ml_tri_area_coords(x1, y1, x2, y2, x3, y3);
        a = ml_sqrt((x2 - x3) * (x2 - x3) + (y2 - y3) * (y2 - y3));
        b = ml_sqrt((x1 - x3) * (x1 - x3) + (y1 - y3) * (y1 - y3));
        c = ml_sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
        R = ml_tri_circumradius(a, b, c);
        r = ml_tri_inradius(a, b, c);
        CHECK(R + 1e-9 * (1.0 + R) >= 2.0 * r, "fx euler inequality R>=2r");
        CHECK(ml_fabs(area - ml_heron(a, b, c)) < 1e-9 * (1.0 + area), "fx area==heron");
        CHECK(ml_euler_line_check(x1, y1, x2, y2, x3, y3, 1e-6) == 1, "fx euler line");
        CHECK(ml_nine_point_check(x1, y1, x2, y2, x3, y3, 1e-6) == 1, "fx nine-point");
        {
            /* Random point set hull: containment + bbox area bound. */
            double xs[12], ys[12];
            int out[12];
            int hn;
            int j;
            int ok = 1;
            double ha;
            for (j = 0; j < 12; j++) {
                xs[j] = 10.0 * fx_unit();
                ys[j] = 10.0 * fx_unit();
            }
            hn = ml_convex_hull(xs, ys, 12, out, 12);
            CHECK(hn >= 3, "fx hull nondegenerate");
            for (j = 0; j < 12 && ok; j++) {
                if (!ml_point_in_polygon(xs[j], ys[j], xs, ys, 12)) {
                    /* containment is w.r.t. the HULL, rebuilt below */
                    ok = ok;
                }
            }
            {
                double hx[12], hy[12];
                for (j = 0; j < hn; j++) {
                    hx[j] = xs[out[j]];
                    hy[j] = ys[out[j]];
                }
                for (j = 0; j < 12; j++) {
                    if (!ml_point_in_polygon(xs[j], ys[j], hx, hy, hn)) {
                        ok = 0;
                        break;
                    }
                }
                ha = ml_polygon_area(hx, hy, hn);
            }
            CHECK(ok, "fx hull contains set");
            CHECK(ha <= 100.0 + 1e-9, "fx hull in bbox");
        }
        CHECK(ml_power_of_point(1.0, 2.0, 1.0, 2.0, 3.0) == -9.0, "fx power center");
        {
            /* Random rectangle area == w*h. */
            double w = 0.5 + 9.5 * fx_unit();
            double h = 0.5 + 9.5 * fx_unit();
            double rxs[4] = {0.0, w, w, 0.0};
            double rys[4] = {0.0, 0.0, h, h};
            double ra = ml_polygon_area(rxs, rys, 4);
            CHECK(ml_fabs(ra - w * h) < 1e-9 * (1.0 + w * h), "fx rect area");
        }
    }
}

/* S6: doctorate spot checks, 400 iters x 20 = 8000. */
static double fx_banach_t(double x) {
    return 0.5 * x + 0.25;
}

void test_fx_doctorate_spot(void) {
    int it;
    printf("--- Doctorate Spot Checks (400 iterations) ---\n");
    for (it = 0; it < 400; it++) {
        {
            double xs = 0;
            CHECK(ml_banach_fixed_point(fx_banach_t, fx_range(0.0, 1.0), 0.5, 1e-12, 10000, &xs, NULL) == ML_SUCCESS, "fx banach st");
            CHECK(ml_fabs(xs - 0.5) < 1e-9, "fx banach value");
        }
        {
            /* Neumann vs direct 2x2 adjugate (diag in (0,2): row-sum contract). */
            double d1 = 0.5 + fx_unit();
            double d2 = 0.5 + fx_unit();
            double o1 = 0.2 * fx_unit();
            double o2 = 0.2 * fx_unit();
            double A[4] = {d1, o1, o2, d2};
            double S[4] = {0, 0, 0, 0};
            double det = d1 * d2 - o1 * o2;
            double res;
            CHECK(ml_neumann_inverse(A, 2, 100, S) == ML_SUCCESS, "fx neumann st");
            res = ml_fabs(A[0] * S[0] + A[1] * S[2] - 1.0) + ml_fabs(A[0] * S[1] + A[1] * S[3]) +
                  ml_fabs(A[2] * S[0] + A[3] * S[2]) + ml_fabs(A[2] * S[1] + A[3] * S[3] - 1.0);
            CHECK(res < 1e-9, "fx neumann residual");
            (void)det;
        }
        {
            /* Hilbert: fixed full-rank V, random y, residual orthogonality. */
            double V[8] = {1, 0, 0, 1, 1, 1, 2, 1};
            double y[4];
            double cc[2] = {0, 0};
            double o1 = 0.0;
            int j;
            for (j = 0; j < 4; j++) y[j] = fx_range(-10.0, 10.0);
            CHECK(ml_hilbert_project(V, y, 4, 2, cc) == ML_SUCCESS, "fx hilbert st");
            for (j = 0; j < 2; j++) {
                double dot = 0.0;
                int k;
                for (k = 0; k < 4; k++) dot += V[k * 2 + j] * (y[k] - (V[k * 2] * cc[0] + V[k * 2 + 1] * cc[1]));
                o1 += (dot < 0.0 ? -dot : dot);
            }
            CHECK(o1 < 1e-9, "fx hilbert orthogonality");
        }
        {
            /* Cycle C_n: betti (1,1), chi 0 (vv holds n + 2n ints). */
            int n = 3 + (int)fx_below(8);
            int vv[30], off[24], dm[23];
            int b[3] = {0, 0, 0};
            int j;
            for (j = 0; j < n; j++) vv[j] = j;
            for (j = 0; j < n - 1; j++) {
                vv[n + 2 * j] = j;
                vv[n + 2 * j + 1] = j + 1;
            }
            /* Closing edge ascending: (0, n-1), not (n-1, 0). */
            vv[n + 2 * (n - 1)] = 0;
            vv[n + 2 * (n - 1) + 1] = n - 1;
            for (j = 0; j <= 2 * n; j++) off[j] = (j <= n) ? j : n + 2 * (j - n);
            for (j = 0; j < n; j++) dm[j] = 0;
            for (j = n; j < 2 * n; j++) dm[j] = 1;
            CHECK(ml_simplicial_betti_f2(n, vv, off, dm, 2 * n, b, 1) == ML_SUCCESS, "fx cycle st");
            CHECK(b[0] == 1 && b[1] == 1, "fx cycle betti");
        }
        {
            /* Bareiss det vs naive double det, tiny ints. */
            long long M[9];
            double nd;
            long long dt = 0;
            int j;
            for (j = 0; j < 9; j++) M[j] = (long long)fx_below(11) - 5;
            nd = (double)M[0] * ((double)M[4] * (double)M[8] - (double)M[5] * (double)M[7]) -
                 (double)M[1] * ((double)M[3] * (double)M[8] - (double)M[5] * (double)M[6]) +
                 (double)M[2] * ((double)M[3] * (double)M[7] - (double)M[4] * (double)M[6]);
            CHECK(ml_bareiss_det(M, 3, &dt) == ML_SUCCESS, "fx bareiss st");
            CHECK((double)dt == nd, "fx bareiss==naive");
        }
        {
            /* Sturm count in range + Descartes dominates it. */
            long long p[4];
            int deg = 1 + (int)fx_below(3);
            int j;
            int cnt;
            int bound;
            for (j = 0; j <= deg; j++) p[j] = (long long)fx_below(7) - 3;
            if (p[deg] == 0) p[deg] = 1;
            cnt = ml_sturm_count(p, deg, -10.0, 10.0);
            bound = ml_descartes_bound(p, deg);
            CHECK(cnt >= 0 && cnt <= deg, "fx sturm range");
            /* Descartes bounds POSITIVE roots only: compare on (0,10).
             * Cauchy: |roots| <= 1+3 = 4 here, so (0,10) catches them all. */
            CHECK(bound >= ml_sturm_count(p, deg, 0.0, 10.0), "fx descartes dominates");
        }
        {
            /* Rodrigues: orthogonality + det 1 + norm preservation. */
            double ax[3];
            double norm;
            double R[9];
            double x[3];
            double y[3];
            int j;
            double o1 = 0.0;
            double det;
            for (j = 0; j < 3; j++) ax[j] = fx_range(-1.0, 1.0);
            norm = ml_sqrt(ax[0] * ax[0] + ax[1] * ax[1] + ax[2] * ax[2]);
            if (norm < 1e-6) {
                ax[0] = 0;
                ax[1] = 0;
                ax[2] = 1;
                norm = 1.0;
            }
            for (j = 0; j < 3; j++) ax[j] /= norm;
            CHECK(ml_rodrigues_exp(ax, fx_range(-3.2, 3.2), R) == ML_SUCCESS, "fx rodrigues st");
            for (j = 0; j < 3; j++) {
                int k;
                for (k = 0; k < 3; k++) {
                    double dot = R[j * 3] * R[k * 3] + R[j * 3 + 1] * R[k * 3 + 1] + R[j * 3 + 2] * R[k * 3 + 2];
                    double want = (j == k) ? 1.0 : 0.0;
                    o1 += (dot > want) ? dot - want : want - dot;
                }
            }
            det = R[0] * (R[4] * R[8] - R[5] * R[7]) - R[1] * (R[3] * R[8] - R[5] * R[6]) + R[2] * (R[3] * R[7] - R[4] * R[6]);
            CHECK(o1 < 1e-12, "fx rodrigues orthogonal");
            CHECK(ml_fabs(det - 1.0) < 1e-12, "fx rodrigues det");
            for (j = 0; j < 3; j++) x[j] = fx_range(-5.0, 5.0);
            for (j = 0; j < 3; j++) y[j] = R[j * 3] * x[0] + R[j * 3 + 1] * x[1] + R[j * 3 + 2] * x[2];
            {
                double nx = ml_sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
                double ny = ml_sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
                CHECK(ml_fabs(nx - ny) < 1e-9 * (1.0 + nx), "fx rotation preserves norm");
            }
        }
        {
            /* Monic cubic from 3 distinct small integer roots. */
            int r0 = (int)fx_below(7) - 3;
            int r1 = r0 + 1 + (int)fx_below(3);
            int r2 = r1 + 1 + (int)fx_below(3);
            double a = 1.0, b = -((double)r0 + (double)r1 + (double)r2);
            double c = (double)r0 * (double)r1 + (double)r0 * (double)r2 + (double)r1 * (double)r2;
            double d = -((double)r0 * (double)r1 * (double)r2);
            double q0 = 0, q1 = 0, q2 = 0;
            double s1 = 0, s2 = 0, s3 = 0;
            int cnt = ml_cubic(a, b, c, d, &q0, &q1, &q2);
            double qs[3] = {q0, q1, q2};
            int j;
            int ok = (cnt == 3);
            CHECK(cnt == 3, "fx cubic finds 3");
            for (j = 0; j < 3; j++) {
                double ev = ((a * qs[j] + b) * qs[j] + c) * qs[j] + d;
                if (ml_fabs(ev) > 1e-6 * (1.0 + ml_fabs(qs[j]) * ml_fabs(qs[j]) * ml_fabs(qs[j]))) ok = 0;
            }
            CHECK(ok, "fx cubic roots verify");
            CHECK(ml_vieta_cubic_sums(a, b, c, d, &s1, &s2, &s3) == ML_SUCCESS, "fx vieta st");
            CHECK(ml_fabs(s1 - (double)(r0 + r1 + r2)) < 1e-9, "fx vieta sums");
        }
    }
}

/* S7: stats/info invariants, 600 iters x 8 = 4800. */
void test_fx_stats_info(void) {
    int it;
    printf("--- Stats/Info Invariants (600 iterations) ---\n");
    for (it = 0; it < 600; it++) {
        double p[4];
        double q[4];
        double s = 0.0;
        double t = 0.0;
        int j;
        double H;
        double KL;
        double CE;
        for (j = 0; j < 4; j++) {
            p[j] = 0.05 + fx_unit();
            q[j] = 0.05 + fx_unit();
            s += p[j];
            t += q[j];
        }
        for (j = 0; j < 4; j++) {
            p[j] /= s;
            q[j] /= t;
        }
        H = ml_entropy(p, 4);
        KL = ml_kl_div(p, q, 4);
        CE = ml_cross_entropy(p, q, 4);
        CHECK(H >= -1e-12, "fx entropy nonneg");
        CHECK(KL >= -1e-9, "fx kl nonneg");
        CHECK(ml_fabs(CE - (H + KL)) < 1e-9 * (1.0 + ml_fabs(H) + ml_fabs(KL)), "fx cross==H+KL");
        {
            double x = fx_range(-10.0, 10.0);
            double L = ml_logistic(x);
            double Lm = ml_logistic(-x);
            double sp = ml_softplus(x);
            CHECK(L > 0.0 && L < 1.0, "fx logistic range");
            CHECK(ml_fabs(L + Lm - 1.0) < 1e-12, "fx logistic symmetry");
            CHECK(sp >= x && sp >= 0.0, "fx softplus dominates");
        }
        {
            /* chi-square roundtrip, moderate regime. */
            int k = 2 + (int)fx_below(4);
            double x = 0.5 + 9.5 * fx_unit();
            double pr = ml_chi2_cdf(x, k);
            double x2 = ml_chi2_inv(pr, k);
            CHECK(ml_fabs(x2 - x) / x < 1e-9, "fx chi2 roundtrip");
        }
        {
            /* Sphere projection lands on the unit sphere. */
            double v[3];
            double w[3] = {0, 0, 0};
            double nrm;
            for (j = 0; j < 3; j++) v[j] = fx_range(-5.0, 5.0);
            if (ml_proj_sphere(v, w, 3) == ML_SUCCESS) {
                nrm = ml_sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
                CHECK(ml_fabs(nrm - 1.0) < 1e-12, "fx sphere unit");
            } else {
                CHECK(0, "fx sphere projects");
            }
        }
    }
}

/* S8: spectral solvers + ODE + cubic, 300 iters x 6 = 1800. */
void test_fx_spectral_ode(void) {
    int it;
    printf("--- Spectral/ODE (300 iterations) ---\n");
    for (it = 0; it < 300; it++) {
        double M[9];
        double A[9];
        double b[3] = {1.0, 2.0, 3.0};
        double x[3] = {0, 0, 0};
        int r;
        int c;
        for (r = 0; r < 9; r++) M[r] = fx_range(-2.0, 2.0);
        for (r = 0; r < 3; r++) {
            for (c = 0; c < 3; c++) {
                double ss = 0.0;
                int k;
                for (k = 0; k < 3; k++) ss += M[k * 3 + r] * M[k * 3 + c];
                A[r * 3 + c] = ss;
            }
        }
        A[0] += 3.0;
        A[4] += 3.0;
        A[8] += 3.0;
        {
            ml_tensor_view_t Av = {A, 3, 3};
            double res = 0.0;
            CHECK(ml_cg_solve(Av, b, x, 100, 1e-12) == ML_SUCCESS, "fx cg st");
            for (r = 0; r < 3; r++) {
                double ss = 0.0;
                for (c = 0; c < 3; c++) ss += A[r * 3 + c] * x[c];
                res += (ss > b[r]) ? ss - b[r] : b[r] - ss;
            }
            CHECK(res < 1e-9, "fx cg residual");
        }
        {
            ml_tensor_view_t Av = {A, 3, 3};
            double re[3] = {0, 0, 0}, im[3] = {0, 0, 0};
            double tr = A[0] + A[4] + A[8];
            double mi = 0.0;
            CHECK(ml_qr_iter_eig(Av, re, im, 500, 1e-10) == ML_SUCCESS, "fx qr st");
            mi = ml_fabs(im[0]) + ml_fabs(im[1]) + ml_fabs(im[2]);
            CHECK(mi < 1e-9, "fx qr real spectrum");
            CHECK(ml_fabs(tr - (re[0] + re[1] + re[2])) < 1e-9, "fx qr trace");
        }
        {
            double y0 = 0.5 + 1.5 * fx_unit();
            double dts[3] = {0.005, 0.01, 0.02};
            double dt = dts[fx_below(3)];
            int steps = (int)(1.0 / dt + 0.5);
            double y = ml_ode_rk4(fx_ode_exp, 0.0, y0, dt, steps);
            double want = y0 * 2.7182818284590451;
            CHECK(ml_fabs(y - want) / want < 1e-6, "fx rk4 exp");
        }
    }
}

/* S9: transform roundtrips + PDE principles, 400 iters x 7 = 2800. */
void test_fx_roundtrips(void) {
    int it;
    printf("--- Transform Roundtrips (400 iterations) ---\n");
    for (it = 0; it < 400; it++) {
        int pws[3] = {8, 16, 32};
        int n = pws[fx_below(3)];
        double x[32];
        int j;
        for (j = 0; j < n; j++) x[j] = fx_range(-5.0, 5.0);
        {
            double avg[1] = {0};
            double det[31] = {0};
            double y[32] = {0};
            double e = 0.0;
            ml_haar_fwt(x, avg, det, n);
            ml_haar_iwt(avg, det, y, n);
            for (j = 0; j < n; j++) {
                double d = (y[j] > x[j]) ? y[j] - x[j] : x[j] - y[j];
                if (d > e) e = d;
            }
            CHECK(e < 1e-12, "fx haar roundtrip");
        }
        {
            double X[32] = {0};
            double z[32] = {0};
            double e = 0.0;
            ml_dct2(x, X, n);
            ml_dct3(X, z, n);
            for (j = 0; j < n; j++) {
                double d = (z[j] > x[j]) ? z[j] - x[j] : x[j] - z[j];
                if (d > e) e = d;
            }
            CHECK(e < 1e-9, "fx dct roundtrip");
        }
        {
            cplx f[32];
            double e = 0.0;
            for (j = 0; j < n; j++) {
                f[j].real = x[j];
                f[j].imag = 0.0;
            }
            ml_fft_execute(f, n);
            ml_ifft_execute(f, n);
            for (j = 0; j < n; j++) {
                double d = ml_fabs(f[j].real - x[j]) + ml_fabs(f[j].imag);
                if (d > e) e = d;
            }
            CHECK(e < 1e-12, "fx fft roundtrip");
        }
        {
            /* Spline interpolates its own nodes. */
            int ns = 5 + (int)fx_below(8);
            double xs[12], ys[12], m2[12];
            double e = 0.0;
            double acc = (double)(1 + (int)fx_below(3));
            int ok = 1;
            for (j = 0; j < ns; j++) {
                acc += 0.5 + fx_unit();
                xs[j] = acc;
                ys[j] = fx_range(-5.0, 5.0);
            }
            if (ml_cubic_spline_natural(xs, ys, ns, m2) != ML_SUCCESS) ok = 0;
            for (j = 0; j < ns && ok; j++) {
                double d = ml_fabs(ml_cubic_spline_eval(xs, ys, m2, ns, xs[j]) - ys[j]);
                if (d > e) e = d;
            }
            CHECK(ok && e < 1e-12, "fx spline nodes");
        }
        {
            /* Poisson residual for -u''=f + heat max principle. */
            int np = 8 + (int)fx_below(16);
            double f[24], u[24];
            double r = 0.0;
            for (j = 0; j < np; j++) f[j] = fx_range(-10.0, 10.0);
            CHECK(ml_poisson_1d(f, u, np, 0.1) == ML_SUCCESS, "fx poisson st");
            for (j = 1; j < np - 1; j++) {
                double lap = (u[j - 1] - 2.0 * u[j] + u[j + 1]) / 0.01;
                double d = ml_fabs(-lap - f[j]);
                if (d > r) r = d;
            }
            CHECK(r < 1e-9, "fx poisson residual");
        }
        {
            double u0[16], u1[16];
            double mx0 = -1e300, mn0 = 1e300;
            double mx1 = -1e300, mn1 = 1e300;
            int ok = 1;
            for (j = 0; j < 16; j++) {
                u0[j] = fx_range(-5.0, 5.0);
                if (u0[j] > mx0) mx0 = u0[j];
                if (u0[j] < mn0) mn0 = u0[j];
            }
            if (ml_heat_explicit(u0, u1, 16, 0.1, 0.004, 1.0) != ML_SUCCESS) ok = 0;
            for (j = 0; j < 16 && ok; j++) {
                if (u1[j] > mx1) mx1 = u1[j];
                if (u1[j] < mn1) mn1 = u1[j];
            }
            CHECK(ok && mx1 <= mx0 + 1e-12 && mn1 >= mn0 - 1e-12, "fx heat max principle");
        }
    }
}

/* S10: integer edge blizzard, 2000 iters x 5 = 10000. */
void test_fx_integer_blizzard(void) {
    int it;
    printf("--- Integer Edge Blizzard (2000 iterations) ---\n");
    for (it = 0; it < 2000; it++) {
        int64_t a;
        int64_t b;
        uint64_t ua;
        uint64_t ub;
        if (fx_below(20) < 3) {
            static const int64_t specials[8] = {0, 1, -1, 2, -2, 1000, -1000, 999983};
            a = specials[fx_below(8)];
            b = specials[fx_below(8)];
        } else {
            a = (int64_t)fx_below(2001) - 1000;
            b = (int64_t)fx_below(2001) - 1000;
        }
        ua = (a < 0) ? (uint64_t)(-(a + 1)) + 1ULL : (uint64_t)a;
        ub = (b < 0) ? (uint64_t)(-(b + 1)) + 1ULL : (uint64_t)b;
        {
            /* egcd agrees with gcd + Bezout holds for small inputs. */
            int64_t x = 0, y = 0;
            int64_t g = ml_egcd(a, b, &x, &y);
            uint64_t want = ml_gcd(ua, ub);
            int ok = (g == (int64_t)want || (want > (uint64_t)INT64_MAX && g == -1));
            if (ok && a != INT64_MIN && b != INT64_MIN && x != 0 && y != 0) {
                long double lhs = (long double)a * (long double)x + (long double)b * (long double)y;
                ok = (lhs == (long double)g);
            }
            CHECK(ok, "fx egcd agrees");
            {
                int64_t g2 = ml_egcd(a, b, NULL, NULL);
                CHECK(g2 == g, "fx egcd null outs");
            }
        }
        {
            /* modinv roundtrip on coprime pairs. */
            uint64_t m = 2 + fx_below(2000);
            uint64_t av = 1 + fx_below(m - 1);
            uint64_t r = ml_modinv(av, m);
            if (ml_gcd(av, m) == 1) {
                CHECK(r != UINT64_MAX && (av * r) % m == 1, "fx modinv roundtrip");
            } else {
                CHECK(r == UINT64_MAX, "fx modinv non-coprime");
            }
        }
        {
            /* crt2 on constructed-compatible pairs. */
            uint64_t m1 = 2 + fx_below(50);
            uint64_t m2 = 2 + fx_below(50);
            uint64_t x0 = fx_below(2000);
            uint64_t a1 = x0 % m1;
            uint64_t a2 = x0 % m2;
            uint64_t x = 0, l = 0;
            ml_status_t st = ml_crt2(a1, m1, a2, m2, &x, &l);
            CHECK(st == ML_SUCCESS && x % m1 == a1 && x % m2 == a2, "fx crt verify");
        }
        {
            /* nCr symmetry, small exact. */
            int n = (int)fx_below(31);
            int k = (n == 0) ? 0 : (int)fx_below((uint64_t)n + 1U);
            CHECK(ml_ncr(n, k) == ml_ncr(n, n - k), "fx ncr symmetry");
        }
    }
}

int main(int argc, char **argv) {
    unsigned int seed;

    if (argc > 1) {
        seed = (unsigned int)strtoul(argv[1], NULL, 10);
    } else {
        const char *env = getenv("MATHLIB_FUZZ_SEED");
        if (env && *env) seed = (unsigned int)strtoul(env, NULL, 10);
        else seed = (unsigned int)time(NULL);
    }

    srand(seed);
    /* Part 2 stream: mixed from the same seed, independent of rand(). */
    fx_state = (uint64_t)seed * 0x9E3779B97F4A7C15ULL + 0xBF58476D1CE4E5B9ULL;
    if (fx_state == 0) fx_state = 0x9E3779B97F4A7C15ULL;
    printf("MATHLIB_FUZZ_SEED=%u\n", seed);
    printf("Fuzz Seed: %u\n", seed);
    printf("=========================================================\n");
    printf("   MATHLIB v12R2: GOD-MODE DYNAMIC FUZZING GAUNTLET\n");
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

    test_fx_contest_identities();
    test_fx_olympiad_nt();
    test_fx_graph_invariants();
    test_fx_modeling_forms();
    test_fx_geometry_predicates();
    test_fx_doctorate_spot();
    test_fx_stats_info();
    test_fx_spectral_ode();
    test_fx_roundtrips();
    test_fx_integer_blizzard();

    printf("\n=========================================================\n");
    printf("GOD-MODE SUMMARY: %d passed, %d failed\n", passed, failed);
    printf("=========================================================\n");

    return failed > 0 ? 1 : 0;
}
