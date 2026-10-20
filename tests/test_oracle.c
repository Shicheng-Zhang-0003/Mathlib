/*
 * MathLib v12A1: Differential Semantic Test (Oracle Validation)
 * MATHLIB_V12A1_ORACLE_EXPANSION
 *
 * Validates: sin, cos, exp, log, gamma, lgamma, pow, atan, and extensions
 * against mpmath ground truth.
 */
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include "ml_core.h"
#include "ml_trig.h"
#include "ml_exp_log.h"
#include "ml_integral.h"

#ifdef MATHLIB_HAS_ORACLE_DATA
#include "oracle_data.h"
#include "oracle_data_ext.h"
#endif
static uint64_t ulp_distance (double a, double b) {
    uint64_t ia, ib;
    memcpy (&ia, &a, sizeof (uint64_t));
    memcpy (&ib, &b, sizeof (uint64_t));
    if (ia >> 63) ia = 0x8000000000000000ULL - ia;
    if (ib >> 63) ib = 0x8000000000000000ULL - ib;
    return ia > ib ? ia - ib : ib - ia;
}
int main () {
#ifdef MATHLIB_HAS_ORACLE_DATA
    int passed = 0;
    int failed = 0;
    int64_t max_ulp = 0;
    double worst_x = 0;
    const char *worst_func = "";

    printf ("=========================================================\n");
    printf ("   MATHLIB v12A1: ORACLE VALIDATION (mpmath ground truth)\n");
    printf ("=========================================================\n");

    /* --- sin --- */
    printf ("--- ml_sin vs mpmath oracle (%d entries) ---\n", oracle_sin_count);
    for (int i = 0; i < oracle_sin_count; i++) {
        double got = ml_sin (oracle_sin [i].input);
        double exp_val = oracle_sin [i].expected;
        if (ml_isnan (exp_val) && ml_isnan (got)) {
            passed++;
            continue;
        }
        if (ml_isinf (exp_val) && ml_isinf (got)) {
            passed++;
            continue;
        }
        int64_t ulp = (int64_t) ulp_distance (got, exp_val);
        if (ulp > max_ulp) {
            max_ulp = ulp;
            worst_x = oracle_sin [i].input;
            worst_func = "sin";
        }
        if (ulp <= 5) {
            passed++;
        } else {
            failed++;
            printf ("  [FAIL] sin(%.17e): got %.17e expected %.17e (%lld ULP)\n",
                    oracle_sin [i].input, got, exp_val, (long long) ulp);
        }
    }

    /* --- cos --- */
    printf ("--- ml_cos vs mpmath oracle (%d entries) ---\n", oracle_cos_count);
    for (int i = 0; i < oracle_cos_count; i++) {
        double got = ml_cos (oracle_cos [i].input);
        double exp_val = oracle_cos [i].expected;
        if (ml_isnan (exp_val) && ml_isnan (got)) {
            passed++;
            continue;
        }
        if (ml_isinf (exp_val) && ml_isinf (got)) {
            passed++;
            continue;
        }
        int64_t ulp = (int64_t) ulp_distance (got, exp_val);
        if (ulp > max_ulp) {
            max_ulp = ulp;
            worst_x = oracle_cos [i].input;
            worst_func = "cos";
        }
        if (ulp <= 5) {
            passed++;
        } else {
            failed++;
            printf ("  [FAIL] cos(%.17e): got %.17e expected %.17e (%lld ULP)\n",
                    oracle_cos [i].input, got, exp_val, (long long) ulp);
        }
    }

    /* --- exp --- */
    printf ("--- ml_exp vs mpmath oracle (%d entries) ---\n", oracle_exp_count);
    for (int i = 0; i < oracle_exp_count; i++) {
        double got = ml_exp (oracle_exp [i].input);
        double exp_val = oracle_exp [i].expected;
        if (ml_isnan (exp_val) && ml_isnan (got)) {
            passed++;
            continue;
        }
        if (ml_isinf (exp_val) && ml_isinf (got)) {
            passed++;
            continue;
        }
        if (exp_val == 0.0 && got == 0.0) {
            passed++;
            continue;
        }
        int64_t ulp = (int64_t) ulp_distance (got, exp_val);
        if (ulp > max_ulp) {
            max_ulp = ulp;
            worst_x = oracle_exp [i].input;
            worst_func = "exp";
        }
        if (ulp <= 5) {
            passed++;
        } else {
            failed++;
            printf ("  [FAIL] exp(%.17e): got %.17e expected %.17e (%lld ULP)\n",
                    oracle_exp [i].input, got, exp_val, (long long) ulp);
        }
    }

    /* --- log --- */
    printf ("--- ml_log vs mpmath oracle (%d entries) ---\n", oracle_log_count);
    for (int i = 0; i < oracle_log_count; i++) {
        double got = ml_log (oracle_log [i].input);
        double exp_val = oracle_log [i].expected;
        if (ml_isnan (exp_val) && ml_isnan (got)) {
            passed++;
            continue;
        }
        if (ml_isinf (exp_val) && ml_isinf (got)) {
            passed++;
            continue;
        }
        int64_t ulp = (int64_t) ulp_distance (got, exp_val);
        if (ulp > max_ulp) {
            max_ulp = ulp;
            worst_x = oracle_log [i].input;
            worst_func = "log";
        }
        if (ulp <= 5) {
            passed++;
        } else {
            failed++;
            printf ("  [FAIL] log(%.17e): got %.17e expected %.17e (%lld ULP)\n",
                    oracle_log [i].input, got, exp_val, (long long) ulp);
        }
    }

    /* --- gamma (v12A1) --- */
    printf ("--- ml_gamma_new vs mpmath oracle (%d entries) ---\n", oracle_gamma_count);
    for (int i = 0; i < oracle_gamma_count; i++) {
        double got = ml_gamma_new (oracle_gamma [i].input);
        double exp_val = oracle_gamma [i].expected;
        if (ml_isnan (exp_val) && ml_isnan (got)) {
            passed++;
            continue;
        }
        if (ml_isinf (exp_val) && ml_isinf (got)) {
            passed++;
            continue;
        }
        if (exp_val == 0.0 && got == 0.0) {
            passed++;
            continue;
        }
        int64_t ulp = (int64_t) ulp_distance (got, exp_val);
        if (ulp > max_ulp) {
            max_ulp = ulp;
            worst_x = oracle_gamma [i].input;
            worst_func = "gamma";
        }
        if (ulp <= 5) {
            passed++;
        } else {
            failed++;
            printf ("  [FAIL] gamma(%.17e): got %.17e expected %.17e (%lld ULP)\n",
                    oracle_gamma [i].input, got, exp_val, (long long) ulp);
        }
    }

    /* --- lgamma (v12A1) --- */
    printf ("--- ml_lgamma vs mpmath oracle (%d entries) ---\n", oracle_lgamma_count);
    for (int i = 0; i < oracle_lgamma_count; i++) {
        double got = ml_lgamma (oracle_lgamma [i].input);
        double exp_val = oracle_lgamma [i].expected;
        if (ml_isnan (exp_val) && ml_isnan (got)) {
            passed++;
            continue;
        }
        if (ml_isinf (exp_val) && ml_isinf (got)) {
            passed++;
            continue;
        }
        int64_t ulp = (int64_t) ulp_distance (got, exp_val);
        if (ulp > max_ulp) {
            max_ulp = ulp;
            worst_x = oracle_lgamma [i].input;
            worst_func = "lgamma";
        }
        if (ulp <= 5) {
            passed++;
        } else {
            failed++;
            printf ("  [FAIL] lgamma(%.17e): got %.17e expected %.17e (%lld ULP)\n",
                    oracle_lgamma [i].input, got, exp_val, (long long) ulp);
        }
    }

    /* --- pow (v12A1) --- */
    printf ("--- ml_pow vs mpmath oracle (%d entries) ---\n", oracle_pow_count);
    for (int i = 0; i < oracle_pow_count; i++) {
        double got = ml_pow (oracle_pow [i].base, oracle_pow [i].exponent);
        double exp_val = oracle_pow [i].expected;
        if (ml_isnan (exp_val) && ml_isnan (got)) {
            passed++;
            continue;
        }
        if (ml_isinf (exp_val) && ml_isinf (got)) {
            passed++;
            continue;
        }
        if (exp_val == 0.0 && got == 0.0) {
            passed++;
            continue;
        }
        int64_t ulp = (int64_t) ulp_distance (got, exp_val);
        if (ulp > max_ulp) {
            max_ulp = ulp;
            worst_x = oracle_pow [i].base;
            worst_func = "pow";
        }
        if (ulp <= 5) {
            passed++;
        } else {
            failed++;
            printf ("  [FAIL] pow(%.17e, %.17e): got %.17e expected %.17e (%lld ULP)\n",
                    oracle_pow [i].base, oracle_pow [i].exponent, got, exp_val, (long long) ulp);
        }
    }

    /* --- R6 ULP-push expansion: 19 families at measured 0 ULP ---
     * Same <=5 ULP gate, same worst-case tracking. Zero/NaN/Inf inputs
     * pass through on exact agreement (mirrors the blocks above). */
#define ORACLE_UNARY(fname, func, arr, count)                                            \
    do {                                                                                 \
        printf ("--- " fname " vs mpmath oracle (%d entries) ---\n", count);             \
        for (int i = 0; i < count; i++) {                                                \
            double got = func ((arr) [i].input);                                         \
            double exp_val = (arr) [i].expected;                                         \
            if (ml_isnan (exp_val) && ml_isnan (got)) {                                  \
                passed++;                                                                \
                continue;                                                                \
            }                                                                            \
            if (ml_isinf (exp_val) && ml_isinf (got)) {                                  \
                passed++;                                                                \
                continue;                                                                \
            }                                                                            \
            if (exp_val == 0.0 && got == 0.0) {                                          \
                passed++;                                                                \
                continue;                                                                \
            }                                                                            \
            {                                                                            \
                int64_t ulp = (int64_t) ulp_distance (got, exp_val);                     \
                if (ulp > max_ulp) {                                                     \
                    max_ulp = ulp;                                                       \
                    worst_x = (arr) [i].input;                                           \
                    worst_func = fname;                                                  \
                }                                                                        \
                if (ulp <= 5) {                                                          \
                    passed++;                                                            \
                } else {                                                                 \
                    failed++;                                                            \
                    printf ("  [FAIL] %s(%.17e): got %.17e expected %.17e (%lld ULP)\n", \
                            fname, (arr) [i].input, got, exp_val, (long long) ulp);      \
                }                                                                        \
            }                                                                            \
        }                                                                                \
    } while (0)

    ORACLE_UNARY ("asin", ml_asin, oracle_asin, oracle_asin_count);
    ORACLE_UNARY ("acos", ml_acos, oracle_acos, oracle_acos_count);
    ORACLE_UNARY ("asinh", ml_asinh, oracle_asinh, oracle_asinh_count);
    ORACLE_UNARY ("acosh", ml_acosh, oracle_acosh, oracle_acosh_count);
    ORACLE_UNARY ("atanh", ml_atanh, oracle_atanh, oracle_atanh_count);
    ORACLE_UNARY ("atan", ml_atan, oracle_atan, oracle_atan_count);
    ORACLE_UNARY ("log10", ml_log10, oracle_log10, oracle_log10_count);
    ORACLE_UNARY ("log2", ml_log2, oracle_log2, oracle_log2_count);
    ORACLE_UNARY ("cbrt", ml_cbrt, oracle_cbrt, oracle_cbrt_count);
    ORACLE_UNARY ("erf", ml_erf, oracle_erf, oracle_erf_count);
    ORACLE_UNARY ("sinh", ml_sinh, oracle_sinh, oracle_sinh_count);
    ORACLE_UNARY ("tanh", ml_tanh, oracle_tanh, oracle_tanh_count);
    ORACLE_UNARY ("sinpi", ml_sinpi, oracle_sinpi, oracle_sinpi_count);
    ORACLE_UNARY ("cospi", ml_cospi, oracle_cospi, oracle_cospi_count);
    ORACLE_UNARY ("tan", ml_tan, oracle_tan, oracle_tan_count);
    ORACLE_UNARY ("exp2", ml_exp2, oracle_exp2, oracle_exp2_count);
    ORACLE_UNARY ("bessel_i0", ml_bessel_i0, oracle_bessel_i0, oracle_bessel_i0_count);
    ORACLE_UNARY ("bessel_i1", ml_bessel_i1, oracle_bessel_i1, oracle_bessel_i1_count);
    ORACLE_UNARY ("bessel_j0", ml_bessel_j0, oracle_bessel_j0, oracle_bessel_j0_count);
    ORACLE_UNARY ("bessel_j1", ml_bessel_j1, oracle_bessel_j1, oracle_bessel_j1_count);

    /* --- Large-argument domain checks (Payne-Hanek) --- */ printf ("--- Large-argument domain checks (Payne-Hanek) ---\n");
    {
        double large_vals [] = {1e10, 1e15, 1e20, 1e50, 1e100, 1e200, 1e300};
        for (int i = 0; i < 7; i++) {
            double s = ml_sin (large_vals [i]);
            double c = ml_cos (large_vals [i]);
            if (!ml_isnan (s) && !ml_isinf (s) && !ml_isnan (c) && !ml_isinf (c)) {
                passed += 2;
            } else {
                failed += 2;
                printf ("  [FAIL] sin/cos(%.1e) not finite\n", large_vals [i]);
            }
            /* Pythagorean identity */
            double pyth = s * s + c * c;
            if (ml_fabs (pyth - 1.0) < 1e-12) {
                passed++;
            } else {
                failed++;
                printf ("  [FAIL] sin^2+cos^2 at %.1e = %.17e\n", large_vals [i], pyth);
            }
        }
    }

    printf ("\n=========================================================\n");
    printf ("ORACLE SUMMARY: %d passed, %d failed\n", passed, failed);
    printf ("Worst case: %s(%.17e) = %lld ULP\n", worst_func, worst_x, (long long) max_ulp);
    if (max_ulp <= 5) {
        printf ("VERIFIED: All functions within <= 5 ULP of mpmath ground truth.\n");
    } else {
        printf ("WARNING: Maximum ULP error exceeds 5. Review worst case above.\n");
    }
    printf ("=========================================================\n");
    return failed > 0 ? 1 : 0;
#else
    printf ("=========================================================\n");
    printf ("   MATHLIB v12A1: ORACLE VALIDATION (SKIPPED)\n");
    printf ("=========================================================\n");
    printf ("Oracle data not found. Generate it first:\n");
    printf ("  pip install mpmath\n");
    printf ("  python3 scripts/oracles/generate_oracles.py\n");
    printf ("Then recompile with: -DMATHLIB_HAS_ORACLE_DATA\n");
    printf ("=========================================================\n");
    return 0;
#endif
}
