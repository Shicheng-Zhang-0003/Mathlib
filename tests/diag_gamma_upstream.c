/* Diagnostic: trace gamma(0.001) error injection point */
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "ml_core.h"
#include "ml_exp_log.h"
#include "ml_integral.h"

static uint64_t ulp_distance(double a, double b) {
    uint64_t ia, ib;
    memcpy(&ia, &a, sizeof(uint64_t));
    memcpy(&ib, &b, sizeof(uint64_t));
    if (ia >> 63) ia = 0x8000000000000000ULL - ia;
    if (ib >> 63) ib = 0x8000000000000000ULL - ib;
    return ia > ib ? ia - ib : ib - ia;
}

int main(void) {
    double x = 1.00000000000000002e-03;  /* 0.001 */
    double oracle_gamma = 9.99423772484595474e+02;
    double oracle_lgamma = 6.90775527898213682e+00;

    printf("=== GAMMA(0.001) UPSTREAM DIAGNOSTIC ===\n\n");

    /* Step 1: lgamma via ml_lgamma */
    double lg = ml_lgamma(x);
    printf("1. ml_lgamma(0.001)       = %.17e\n", lg);
    printf("   oracle lgamma          = %.17e\n", oracle_lgamma);
    printf("   ULP distance           = %llu\n\n",
           (unsigned long long)ulp_distance(lg, oracle_lgamma));

    /* Step 2: gamma via ml_gamma_new (the failing path) */
    double g = ml_gamma_new(x);
    printf("2. ml_gamma_new(0.001)    = %.17e\n", g);
    printf("   oracle gamma           = %.17e\n", oracle_gamma);
    printf("   ULP distance           = %llu\n\n",
           (unsigned long long)ulp_distance(g, oracle_gamma));

    /* Step 3: exp(ml_lgamma) using LIBC exp (bypasses ml_exp entirely) */
    double g_libc = exp(lg);
    printf("3. libc exp(ml_lgamma)    = %.17e\n", g_libc);
    printf("   oracle gamma           = %.17e\n", oracle_gamma);
    printf("   ULP distance           = %llu\n\n",
           (unsigned long long)ulp_distance(g_libc, oracle_gamma));

    /* Step 4: exp(ml_lgamma) using ml_exp (the suspect) */
    double g_ml = ml_exp(lg);
    printf("4. ml_exp(ml_lgamma)      = %.17e\n", g_ml);
    printf("   oracle gamma           = %.17e\n", oracle_gamma);
    printf("   ULP distance           = %llu\n\n",
           (unsigned long long)ulp_distance(g_ml, oracle_gamma));

    /* Step 5: trace ml_exp internals */
    double ln2 = 0.693147180559945309417;
    double n_d = lg / ln2;
    double n_r = ml_round(n_d);
    int n = (int)n_r;
    printf("5. ml_exp internals for input %.17e:\n", lg);
    printf("   lg / ln2               = %.17e\n", n_d);
    printf("   ml_round(lg / ln2)     = %.1f  (n = %d)\n", n_r, n);

    /* Cody-Waite reduction */
    double LN2_HI = 6.93147180369123816490e-01;
    double LN2_LO = 1.90821492927058500170e-10;
    double r = fma(-n_r, LN2_HI, lg);
    r = fma(-n_r, LN2_LO, r);
    printf("   reduced arg r          = %.17e\n", r);
    printf("   |r|                    = %.17e\n\n", fabs(r));

    /* Step 6: compare ml_exp vs libc exp at the SAME input */
    double test_inputs[] = {6.90775527898213682, 6.9, 7.0, 6.0, 5.0, 10.0};
    printf("6. ml_exp vs libc exp at various inputs:\n");
    for (int i = 0; i < 6; i++) {
        double ti = test_inputs[i];
        double me = ml_exp(ti);
        double le = exp(ti);
        printf("   x=%.17e  ml_exp=%.17e  libc=%.17e  ULP=%llu\n",
               ti, me, le,
               (unsigned long long)ulp_distance(me, le));
    }

    printf("\n=== DIAGNOSIS ===\n");
    unsigned long long ulp_libc = (unsigned long long)ulp_distance(g_libc, oracle_gamma);
    unsigned long long ulp_ml   = (unsigned long long)ulp_distance(g_ml, oracle_gamma);
    if (ulp_libc <= 5 && ulp_ml > 5) {
        printf("ERROR IS IN ml_exp (or ml_ldexp_pure).\n");
        printf("libc exp matches oracle, ml_exp does not.\n");
    } else if (ulp_libc > 5 && ulp_ml > 5) {
        printf("ERROR IS IN THE LGAMMA DD CHAIN.\n");
        printf("Both libc exp and ml_exp disagree with oracle.\n");
        printf("The lgamma value fed to exp is wrong.\n");
    } else {
        printf("BOTH MATCH. Error may be in ml_exp_dd FMA correction.\n");
    }
    printf("=================\n");
    return 0;
}
