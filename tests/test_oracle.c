/*
 * MathLib v11S: Differential Semantic Test (Oracle Validation)
 * Compares ml_sin, ml_cos, ml_exp, ml_log against mpmath ground truth.
 *
 * Prerequisite: Run `python3 scripts/oracles/generate_oracles.py` first
 * to generate tests/oracle_data.h
 */
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include "ml_core.h"
#include "ml_trig.h"
#include "ml_exp_log.h"

#ifdef MATHLIB_HAS_ORACLE_DATA
#include "oracle_data.h"
#endif

/* Exact ULP distance calculator (handles sign-bit crossing) */
static uint64_t ulp_distance(double a, double b) {
    uint64_t ia, ib;
    memcpy(&ia, &a, sizeof(uint64_t));
    memcpy(&ib, &b, sizeof(uint64_t));
    if (ia >> 63) ia = 0x8000000000000000ULL - ia;
    if (ib >> 63) ib = 0x8000000000000000ULL - ib;
    return ia > ib ? ia - ib : ib - ia;
}

int main() {
#ifdef MATHLIB_HAS_ORACLE_DATA
    int passed = 0;
    int failed = 0;
    int64_t max_ulp = 0;
    double worst_x = 0;
    const char* worst_func = "";

    printf("=========================================================\n");
    printf("   MATHLIB v11S: ORACLE VALIDATION (mpmath ground truth)\n");
    printf("=========================================================\n");

    /* Validate sin */
    printf("--- ml_sin vs mpmath oracle (%d entries) ---\n", oracle_sin_count);
    for (int i = 0; i < oracle_sin_count; i++) {
        double got = ml_sin(oracle_sin[i].input);
        double exp_val = oracle_sin[i].expected;
        if (ml_isnan(exp_val) && ml_isnan(got)) { passed++; continue; }
        if (ml_isinf(exp_val) && ml_isinf(got)) { passed++; continue; }
        int64_t ulp = ulp_distance(got, exp_val);
        if (ulp > max_ulp) { max_ulp = ulp; worst_x = oracle_sin[i].input; worst_func = "sin"; }
        if (ulp <= 5) { passed++; }
        else {
            failed++;
            printf("  [FAIL] sin(%.17e): got %.17e expected %.17e (%lld ULP)\n",
                oracle_sin[i].input, got, exp_val, (long long)ulp);
        }
    }

    /* Validate cos */
    printf("--- ml_cos vs mpmath oracle (%d entries) ---\n", oracle_cos_count);
    for (int i = 0; i < oracle_cos_count; i++) {
        double got = ml_cos(oracle_cos[i].input);
        double exp_val = oracle_cos[i].expected;
        if (ml_isnan(exp_val) && ml_isnan(got)) { passed++; continue; }
        if (ml_isinf(exp_val) && ml_isinf(got)) { passed++; continue; }
        int64_t ulp = ulp_distance(got, exp_val);
        if (ulp > max_ulp) { max_ulp = ulp; worst_x = oracle_cos[i].input; worst_func = "cos"; }
        if (ulp <= 5) { passed++; }
        else {
            failed++;
            printf("  [FAIL] cos(%.17e): got %.17e expected %.17e (%lld ULP)\n",
                oracle_cos[i].input, got, exp_val, (long long)ulp);
        }
    }

    /* Validate exp */
    printf("--- ml_exp vs mpmath oracle (%d entries) ---\n", oracle_exp_count);
    for (int i = 0; i < oracle_exp_count; i++) {
        double got = ml_exp(oracle_exp[i].input);
        double exp_val = oracle_exp[i].expected;
        if (ml_isnan(exp_val) && ml_isnan(got)) { passed++; continue; }
        if (ml_isinf(exp_val) && ml_isinf(got)) { passed++; continue; }
        if (exp_val == 0.0 && got == 0.0) { passed++; continue; }
        int64_t ulp = ulp_distance(got, exp_val);
        if (ulp > max_ulp) { max_ulp = ulp; worst_x = oracle_exp[i].input; worst_func = "exp"; }
        if (ulp <= 5) { passed++; }
        else {
            failed++;
            printf("  [FAIL] exp(%.17e): got %.17e expected %.17e (%lld ULP)\n",
                oracle_exp[i].input, got, exp_val, (long long)ulp);
        }
    }

    /* Validate log */
    printf("--- ml_log vs mpmath oracle (%d entries) ---\n", oracle_log_count);
    for (int i = 0; i < oracle_log_count; i++) {
        double got = ml_log(oracle_log[i].input);
        double exp_val = oracle_log[i].expected;
        if (ml_isnan(exp_val) && ml_isnan(got)) { passed++; continue; }
        if (ml_isinf(exp_val) && ml_isinf(got)) { passed++; continue; }
        int64_t ulp = ulp_distance(got, exp_val);
        if (ulp > max_ulp) { max_ulp = ulp; worst_x = oracle_log[i].input; worst_func = "log"; }
        if (ulp <= 5) { passed++; }
        else {
            failed++;
            printf("  [FAIL] log(%.17e): got %.17e expected %.17e (%lld ULP)\n",
                oracle_log[i].input, got, exp_val, (long long)ulp);
        }
    }

    
    /* Validate v11S Domain Boundaries (2-term Cody-Waite limits) */
    printf("--- Domain Boundary Checks ---\n");
    double massive_inputs[] = {1e10, 1e15, -1e10, -1e15};
    for (int i = 0; i < 4; i++) {
        double s = ml_sin(massive_inputs[i]);
        double c = ml_cos(massive_inputs[i]);
        // Beyond 1e6, 2-term Cody-Waite loses precision, but it MUST NOT crash or return Inf/NaN
        if (!ml_isnan(s) && !ml_isinf(s) && !ml_isnan(c) && !ml_isinf(c)) { passed += 2; }
        else { failed += 2; printf("  [FAIL] Domain boundary crashed: x=%.1e\n", massive_inputs[i]); }
    }

    // Beyond 1e18, the library MUST fail loudly with NaN to prevent long long UB
    if (ml_isnan(ml_sin(1e50)) && ml_isnan(ml_cos(1e50))) { passed += 2; }
    else { failed += 2; printf("  [FAIL] 1e50 did not safely return NaN\n"); }

    printf("\n=========================================================\n");
    printf("ORACLE SUMMARY: %d passed, %d failed\n", passed, failed);
    printf("Worst case: %s(%.17e) = %lld ULP\n", worst_func, worst_x, (long long)max_ulp);
    if (max_ulp <= 5) {
        printf("✅ VERIFIED: All functions within <= 5 ULP of mpmath ground truth.\n");
    } else {
        printf("⚠️  WARNING: Maximum ULP error exceeds 5. Review worst case above.\n");
    }
    printf("=========================================================\n");
    return failed > 0 ? 1 : 0;

#else
    printf("=========================================================\n");
    printf("   MATHLIB v11S: ORACLE VALIDATION (SKIPPED)\n");
    printf("=========================================================\n");
    printf("Oracle data not found. Generate it first:\n");
    printf("  pip install mpmath\n");
    printf("  python3 scripts/oracles/generate_oracles.py\n");
    printf("Then recompile with: -DMATHLIB_HAS_ORACLE_DATA\n");
    printf("=========================================================\n");
    return 0;
#endif
}
