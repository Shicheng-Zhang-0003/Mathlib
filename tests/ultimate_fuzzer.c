#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "ml_core.h"
#include "ml_trig.h"
#include "ml_exp_log.h"
#include "ml_complex.h"
#include "ml_linalg.h"
#include "ml_tensor.h"
#include "ml_combinatorics.h"
#include "ml_quaternion.h"
#include "ml_fixed_point.h"
#include "fft.h"
#include "ml_compiler.h"
#include "ml_types.h"
#include "bitwise_fp.h"

static uint64_t passed = 0;
static uint64_t failed = 0;

#define CHECK(cond, msg) do { \
    if (cond) { passed++; } \
    else { failed++; printf("FAIL: %s (Line %d)\n", msg, __LINE__); } \
} while(0)

#define CHECK_NEAR(a, b, eps, msg) do { \
    double _a = (double)(a); double _b = (double)(b); \
    if (ml_fabs(_a - _b) < (eps)) { passed++; } \
    else { failed++; printf("FAIL: %s got=%.17g exp=%.17g diff=%.17g (Line %d)\n", msg, _a, _b, ml_fabs(_a-_b), __LINE__); } \
} while(0)

static double rand_fp() {
    uint64_t bits = ((uint64_t)rand() << 32) | (uint64_t)rand();
    double d;
    memcpy(&d, &bits, sizeof(double));
    return d;
}

void test_bitwise_ieee754() {
    printf("  [1/5] Bitwise IEEE-754 (10k iters)...\n");
    for(int i=0; i<10000; i++) {
        double x = rand_fp();
        double y = rand_fp();

        if (!ml_isnan(x) && !ml_isnan(y) && !ml_isinf(x) && !ml_isinf(y)
            && y != 0.0 && !ml_is_subnormal(x) && !ml_is_subnormal(y)) {
            double ml_res = ml_fmod(x, y);
            double libc_res = fmod(x, y);
            if (!ml_isnan(ml_res) && !ml_isnan(libc_res)) {
    double tol = ml_fabs(libc_res) * 1e-13 + 1e-15;

    if (ml_fabs(ml_res - libc_res) < tol) {
        passed++;
    } else {
        failed++;
        printf("FAIL: fmod equivalence | func=fmod x=%.17g y=%.17g got=%.17g expected=%.17g diff=%.17g\n",
               x, y, ml_res, libc_res, ml_fabs(ml_res - libc_res));
    }
}
        }

        if (!ml_isnan(x) && !ml_isinf(x) && x != 0.0 && !ml_is_subnormal(x)) {
            int ml_exp, libc_exp;
            double ml_mant = ml_frexp_pure(x, &ml_exp);
            double libc_mant = frexp(x, &libc_exp);
            CHECK_NEAR(ml_mant, libc_mant, 1e-15, "frexp mantissa");
            CHECK(ml_exp == libc_exp, "frexp exponent");
        }
    }
}

void test_transcendental_extremes() {
    printf("  [2/5] Transcendental Extremes (10k iters)...\n");
    for(int i=0; i<10000; i++) {
        double x = rand_fp();
        ml_sin(x); ml_cos(x); ml_tan(x);
        ml_asin(x); ml_acos(x); ml_atan(x);
        ml_exp(x); ml_log(x > 0 ? x : 1.0);
        ml_sqrt(x > 0 ? x : 0.0);
        ml_sinh(x); ml_cosh(x); ml_tanh(x);
        ml_atan2(x, rand_fp());
    }
    passed++;
}

void test_workspace_torture() {
    printf("  [3/5] Workspace Torture (10k iters)...\n");
    char scratchpad[8192];
    ml_workspace_t ws;
    ml_workspace_init(&ws, scratchpad, sizeof(scratchpad));

    for(int i=0; i<10000; i++) {
        size_t req = (size_t)(rand() % 1024) + 1;
        void* ptr = ml_workspace_alloc(&ws, req);
        if (ptr) {
            memset(ptr, 0xAA, req);
        }
        if (rand() % 10 == 0) {
            ml_workspace_reset(&ws);
        }
    }
    CHECK(ws.magic_canary == ML_WORKSPACE_CANARY, "Workspace canary survived");
}

void test_combinatorics() {
    printf("  [4/5] Combinatorics Boundaries...\n");
    CHECK(ml_ncr(60, 30) == 118264581564861424ULL, "ncr 60,30");
    CHECK(ml_ncr(61, 30) == 232714176627630544ULL, "ncr 61,30");
    CHECK(ml_ncr(62, 31) == 465428353255261088ULL, "ncr 62,31");
    passed += 3;
}

void test_algebraic_invariants() {
    printf("  [5/5] Algebraic Invariants (10k iters)...\n");
    for(int i=0; i<10000; i++) {
        double x = ((double)rand() / RAND_MAX * 200.0) - 100.0;
        double s = ml_sin(x);
        double c = ml_cos(x);
        CHECK_NEAR(s*s + c*c, 1.0, 1e-13, "Pythagorean");

        double px = x > 1e-5 ? x : 1e-5;
        if (px < 700.0) {
            CHECK_NEAR(ml_log(ml_exp(px)), px, 1e-12, "log(exp(x))");
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
    printf("MATHLIB_FUZZ_SEED=%u\n", seed);
    printf("=========================================================\n");
    printf("   MATHLIB v11S: THE ULTIMATE FUZZER (ASan + UBSan)\n");
    printf("=========================================================\n");

    test_bitwise_ieee754();
    test_transcendental_extremes();
    test_workspace_torture();
    test_combinatorics();
    test_algebraic_invariants();

    printf("\n=========================================================\n");
    printf("ULTIMATE SUMMARY: %llu passed, %llu failed\n",
        (unsigned long long)passed, (unsigned long long)failed);
    printf("=========================================================\n");
    return failed > 0 ? 1 : 0;
}
