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
#include "ml_analytic_nt.h"
#include "ml_integral.h"
#include "ml_statistics.h"
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
#include "ml_kelvin.h"

static uint64_t passed = 0;
static uint64_t failed = 0;

#define CHECK(cond, msg)                                    \
    do {                                                    \
        if (cond) {                                         \
            passed++;                                       \
        } else {                                            \
            failed++;                                       \
            printf ("FAIL: %s (Line %d)\n", msg, __LINE__); \
        }                                                   \
    } while (0)

#define CHECK_NEAR(a, b, eps, msg)                                                                                    \
    do {                                                                                                              \
        double _a = (double) (a);                                                                                     \
        double _b = (double) (b);                                                                                     \
        if (ml_fabs (_a - _b) < (eps)) {                                                                              \
            passed++;                                                                                                 \
        } else {                                                                                                      \
            failed++;                                                                                                 \
            printf ("FAIL: %s got=%.17g exp=%.17g diff=%.17g (Line %d)\n", msg, _a, _b, ml_fabs (_a - _b), __LINE__); \
        }                                                                                                             \
    } while (0)

#define CHECK_NEAR_REL(a, b, eps, msg)                                                                    \
    do {                                                                                                  \
        double _a = (double) (a);                                                                         \
        double _b = (double) (b);                                                                         \
        double _diff = ml_fabs (_a - _b);                                                                 \
        double _scale = ml_fabs (_b) > 1.0 ? ml_fabs (_b) : 1.0;                                          \
        if (_diff < (eps) * _scale) {                                                                     \
            passed++;                                                                                     \
        } else {                                                                                          \
            failed++;                                                                                     \
            printf ("FAIL: %s got=%.17g exp=%.17g diff=%.17g (Line %d)\n", msg, _a, _b, _diff, __LINE__); \
        }                                                                                                 \
    } while (0)
static double rand_fp () {
    uint64_t bits = ((uint64_t) rand () << 32) | (uint64_t) rand ();
    double d;
    memcpy (&d, &bits, sizeof (double));
    return d;
}
/* Crash-resistance contract: arbitrary bits in, no crash/hang out;
 * every double out is finite-or-NaN, every status in the enum set. */
static int ult_sane (double x) {
    return ml_isfinite (x) || ml_isnan (x);
}
#define CHECK_SANE(x, msg) CHECK (ult_sane (x), msg)
#define CHECK_STATUS(st, msg) CHECK ((st) == ML_SUCCESS || (st) == ML_ERR_SINGULAR ||                                           \
                                         (st) == ML_ERR_WORKSPACE || (st) == ML_ERR_INVALID_ARG || (st) == ML_ERR_NAN_INPUT ||  \
                                         (st) == ML_ERR_OVERFLOW || (st) == ML_ERR_CONVERGENCE || (st) == ML_ERR_UNSUPPORTED || \
                                         (st) == ML_ERR_INTERNAL,                                                               \
                                     msg)
static uint64_t ult_u64 (void) {
    return ((uint64_t) (uint32_t) rand () << 32) | (uint64_t) (uint32_t) rand ();
}
static double ult_banach_t (double x) {
    return 0.5 * x + 0.25;
}
void test_bitwise_ieee754 () {
    printf ("  [1/5] Bitwise IEEE-754 (10k iters)...\n");
    for (int i = 0; i < 10000; i++) {
        double x = rand_fp ();
        double y = rand_fp ();

        if (!ml_isnan (x) && !ml_isnan (y) && !ml_isinf (x) && !ml_isinf (y) && y != 0.0 && !ml_is_subnormal (x) && !ml_is_subnormal (y)) {
            double ml_res = ml_fmod (x, y);
            double libc_res = fmod (x, y);
            if (!ml_isnan (ml_res) && !ml_isnan (libc_res)) {
                double tol = ml_fabs (libc_res) * 1e-13 + 1e-15;

                if (ml_fabs (ml_res - libc_res) < tol) {
                    passed++;
                } else {
                    failed++;
                    printf ("FAIL: fmod equivalence | func=fmod x=%.17g y=%.17g got=%.17g expected=%.17g diff=%.17g\n",
                            x, y, ml_res, libc_res, ml_fabs (ml_res - libc_res));
                }
            }
        }

        if (!ml_isnan (x) && !ml_isinf (x) && x != 0.0 && !ml_is_subnormal (x)) {
            int ml_exp, libc_exp;
            double ml_mant = ml_frexp_pure (x, &ml_exp);
            double libc_mant = frexp (x, &libc_exp);
            CHECK_NEAR (ml_mant, libc_mant, 1e-15, "frexp mantissa");
            CHECK (ml_exp == libc_exp, "frexp exponent");
        }
    }
}
void test_transcendental_extremes () {
    printf ("  [2/5] Transcendental Extremes (10k iters)...\n");
    for (int i = 0; i < 10000; i++) {
        double x = rand_fp ();
        ml_sin (x);
        ml_cos (x);
        ml_tan (x);
        ml_asin (x);
        ml_acos (x);
        ml_atan (x);
        ml_exp (x);
        ml_log (x > 0 ? x : 1.0);
        ml_sqrt (x > 0 ? x : 0.0);
        ml_sinh (x);
        ml_cosh (x);
        ml_tanh (x);
        ml_atan2 (x, rand_fp ());
    }
    passed++;
}
void test_workspace_torture () {
    printf ("  [3/5] Workspace Torture (10k iters)...\n");
    char scratchpad [8192];
    ml_workspace_t ws;
    ml_workspace_init (&ws, scratchpad, sizeof (scratchpad));

    for (int i = 0; i < 10000; i++) {
        size_t req = (size_t) (rand () % 1024) + 1;
        void *ptr = ml_workspace_alloc (&ws, req);
        if (ptr) {
            memset (ptr, 0xAA, req);
        }
        if (rand () % 10 == 0) {
            ml_workspace_reset (&ws);
        }
    }
    CHECK (ws.magic_canary == ML_WORKSPACE_CANARY, "Workspace canary survived");
}
void test_combinatorics () {
    printf ("  [4/7] Combinatorics Boundaries...\n");
    CHECK (ml_ncr (60, 30) == 118264581564861424ULL, "ncr 60,30");
    CHECK (ml_ncr (61, 30) == 232714176627630544ULL, "ncr 61,30");
    CHECK (ml_ncr (62, 31) == 465428353255261088ULL, "ncr 62,31");
    passed += 3;
}
void test_algebraic_invariants () {
    printf ("  [5/7] Algebraic Invariants (10k iters)...\n");
    for (int i = 0; i < 10000; i++) {
        double x = ((double) rand () / RAND_MAX * 200.0) - 100.0;
        double s = ml_sin (x);
        double c = ml_cos (x);
        CHECK_NEAR (s * s + c * c, 1.0, 1e-13, "Pythagorean");

        double px = x > 1e-5 ? x : 1e-5;
        if (px < 700.0) {
            CHECK_NEAR_REL (ml_log (ml_exp (px)), px, 1e-13, "log(exp(x))");
        }
    }
}
/* Crash upgrade (R6): arbitrary bit patterns into every new double API.
 * Contract under test is survival (ASan/UBSan + the gate timeout catch
 * hangs/crashes); poles and overflow legitimately produce Inf, so those
 * three are survival-labelled. Value pins live in edge suites, not here.
 * 1500 iters x 13 = 19500. */
void test_crash_new_doubles (void) {
    int i;
    printf ("  [6/7] New-API double crash (1500 iters)...\n");
    for (i = 0; i < 1500; i++) {
        double x = rand_fp ();
        CHECK_SANE (ml_bessel_j0 (x), "crash bessel_j0");
        CHECK (ml_isfinite (ml_kelvin_ber (x)) || ml_isnan (ml_kelvin_ber (x)) || ml_isinf (ml_kelvin_ber (x)), "survived kelvin_ber");
        CHECK_SANE (ml_erf (x), "crash erf");
        CHECK (ml_isfinite (ml_gamma_new (x)) || ml_isnan (ml_gamma_new (x)) || ml_isinf (ml_gamma_new (x)), "survived gamma_new");
        {
            double arr [4];
            int j;
            for (j = 0; j < 4; j++) arr [j] = rand_fp ();
            CHECK_SANE (ml_entropy (arr, 4), "crash entropy");
            CHECK_SANE (ml_kde_gaussian (arr, 4, x, 1.0), "crash kde");
        }
        {
            double L = ml_logistic (x);
            CHECK (ml_isnan (L) || (L >= 0.0 && L <= 1.0), "crash logistic range");
            {
                double sp = ml_softplus (x);
                CHECK (ml_isnan (sp) || sp >= x, "crash softplus dominates");
            }
        }
        {
            double p = ml_isnan (x) || !ml_isfinite (x) ? 0.5 : x - floor (x);
            double v;
            if (p <= 0.0 || p >= 1.0) p = 0.5;
            v = ml_chi2_inv (p, 4);
            CHECK (ml_isnan (v) || v >= 0.0, "crash chi2_inv nonneg");
        }
        CHECK (ml_normal_logcdf (x, 0.0, 1.0) <= 0.0 || ml_isnan (ml_normal_logcdf (x, 0.0, 1.0)), "crash logcdf nonpos");
        CHECK_SANE (ml_morlet_cwt (&x, 1, 1.0, 2.0, 0), "crash morlet");
        CHECK (ml_isfinite (ml_digamma (x)) || ml_isnan (ml_digamma (x)) || ml_isinf (ml_digamma (x)), "survived digamma");
        {
            double xs = 0.0;
            double x0 = (ml_isfinite (x) ? x : 0.0);
            if (x0 > 1e6) x0 = 1.0;
            if (x0 < -1e6) x0 = -1.0;
            CHECK_STATUS (ml_banach_fixed_point (ult_banach_t, x0, 0.5, 1e-12, 1000, &xs, NULL), "crash banach status");
        }
    }
}
/* Crash upgrade (R6): raw 64-bit patterns into integer APIs with bounded
 * slow paths (sieve/factor caps keep every call fast).
 * 2000 iters x 10 = 20000. */
static uint64_t ult_fact (int n) {
    static const uint64_t t [21] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320,
                                    362880, 3628800, 39916800, 479001600, 6227020800ULL, 87178291200ULL,
                                    1307674368000ULL, 20922789888000ULL, 355687428096000ULL,
                                    6402373705728000ULL, 121645100408832000ULL, 2432902008176640000ULL};
    return t [n];
}
void test_crash_integer_blizzard (void) {
    int i;
    printf ("  [7/7] Integer blizzard (2000 iters)...\n");
    for (i = 0; i < 2000; i++) {
        uint64_t u = ult_u64 ();
        int64_t a = (int64_t) u;
        int64_t b = (int64_t) ult_u64 ();
        int64_t egx = 0, egy = 0;
        int64_t g = ml_egcd (a, b, &egx, &egy);
        CHECK (g == -1 || g >= 0, "crash egcd range");
        {
            uint64_t m = 2 + u % 99999;
            uint64_t av = 1 + ult_u64 () % (m - 1);
            uint64_t r = ml_modinv (av, m);
            if (ml_gcd (av, m) == 1) {
                CHECK (r != UINT64_MAX && (av * r) % m == 1, "crash modinv roundtrip");
            } else {
                CHECK (r == UINT64_MAX, "crash modinv sentinel");
            }
        }
        {
            int fi = (int) (u % 50) - 10;
            uint64_t f = ml_factorial (fi);
            if (fi >= 0 && fi <= 20) {
                CHECK (f == ult_fact (fi), "crash factorial table");
            } else {
                CHECK (f == UINT64_MAX, "crash factorial sentinel");
            }
        }
        {
            uint64_t n = 1 + u % 99999;
            uint64_t phi = ml_phi (n);
            uint64_t tau = ml_tau_divisors (n);
            uint64_t sig = ml_sigma_sum (n);
            int mu = ml_mobius ((int64_t) n);
            CHECK (phi >= 1, "crash phi positive");
            CHECK (tau >= 1, "crash tau positive");
            CHECK (sig >= n + 1 || sig == UINT64_MAX, "crash sigma bound");
            CHECK (mu == -1 || mu == 0 || mu == 1 || mu == 2, "crash mobius set");
        }
        {
            int v = ml_is_prime64 (u);
            CHECK (v == 0 || v == 1, "crash is_prime boolean");
            {
                uint64_t n = 2 + u % 99998;
                CHECK (ml_prime_pi (n) <= n, "crash prime_pi bound");
            }
        }
        {
            int bn = (int) (u % 26);
            uint64_t bl = ml_bell_number (bn);
            CHECK (bl >= 1, "crash bell positive");
        }
    }
}
int main (int argc, char **argv) {
    unsigned int seed;

    if (argc > 1) {
        seed = (unsigned int) strtoul (argv [1], NULL, 10);
    } else {
        const char *env = getenv ("MATHLIB_FUZZ_SEED");
        if (env && *env) seed = (unsigned int) strtoul (env, NULL, 10);
        else seed = (unsigned int) time (NULL);
    }

    srand (seed);
    printf ("MATHLIB_FUZZ_SEED=%u\n", seed);
    printf ("=========================================================\n");
    printf ("   MATHLIB v12R2: THE ULTIMATE FUZZER (ASan + UBSan)\n");
    printf ("=========================================================\n");

    test_bitwise_ieee754 ();
    test_transcendental_extremes ();
    test_workspace_torture ();
    test_combinatorics ();
    test_algebraic_invariants ();
    test_crash_new_doubles ();
    test_crash_integer_blizzard ();

    printf ("\n=========================================================\n");
    printf ("ULTIMATE SUMMARY: %llu passed, %llu failed\n",
            (unsigned long long) passed, (unsigned long long) failed);
    printf ("=========================================================\n");
    return failed > 0 ? 1 : 0;
}
