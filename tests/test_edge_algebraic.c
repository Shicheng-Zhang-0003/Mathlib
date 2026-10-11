/* Doctorate algebraic pins: resultant, discriminant, Sturm, Descartes,
 * hook, class sizes, Rodrigues, Ei/li/psi. Truth: sympy + closed forms. */
#include "test_harness.h"
#include "ml_algebraic.h"
#include "ml_analytic_nt.h"
#include "ml_core.h"

int main(void) {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Edge Algebraic");

    /* Res(x^2-2, x^3-3) = 1; Res(x^2-1, x^2-3x+2) = 0 (share x=1). */
    {
        long long r = 0;
        long long p[] = {-2, 0, 1}, q[] = {-3, 0, 0, 1};
        ASSERT_TRUE(&ctx, ml_resultant(p, 2, q, 3, &r) == ML_SUCCESS && r == 1, "resultant 1");
        {
            long long a[] = {-1, 0, 1}, b[] = {2, -3, 1};
            ASSERT_TRUE(&ctx, ml_resultant(a, 2, b, 2, &r) == ML_SUCCESS && r == 0, "resultant common root");
        }
        ASSERT_TRUE(&ctx, ml_resultant(NULL, 2, q, 3, &r) == ML_ERR_INVALID_ARG, "resultant null");
    }
    /* disc((x-1)(x-2)(x-3)) = 4; disc(x^2+1) = -4; linear = 1. */
    {
        long long d = 0;
        long long p[] = {-6, 11, -6, 1};
        long long q[] = {1, 0, 1};
        long long l[] = {3, 2};
        ASSERT_TRUE(&ctx, ml_discriminant(p, 3, &d) == ML_SUCCESS && d == 4, "disc cubic");
        ASSERT_TRUE(&ctx, ml_discriminant(q, 2, &d) == ML_SUCCESS && d == -4, "disc x^2+1");
        ASSERT_TRUE(&ctx, ml_discriminant(l, 1, &d) == ML_SUCCESS && d == 1, "disc linear");
    }
    /* Sturm: 3 roots in (0,4); x^2+1 none; x^2-2 one in (0,2). */
    {
        long long p[] = {-6, 11, -6, 1};
        long long q[] = {1, 0, 1};
        long long r[] = {-2, 0, 1};
        ASSERT_TRUE(&ctx, ml_sturm_count(p, 3, 0.0, 4.0) == 3, "sturm cubic");
        ASSERT_TRUE(&ctx, ml_sturm_count(p, 3, 1.5, 2.5) == 1, "sturm isolate [2]");
        ASSERT_TRUE(&ctx, ml_sturm_count(q, 2, -10.0, 10.0) == 0, "sturm none");
        ASSERT_TRUE(&ctx, ml_sturm_count(r, 2, 0.0, 2.0) == 1, "sturm sqrt2");
        ASSERT_TRUE(&ctx, ml_sturm_count(p, 3, 4.0, 0.0) == -1, "sturm flipped");
    }
    ASSERT_TRUE(&ctx, ml_descartes_bound((long long[]){-6, 11, -6, 1}, 3) == 3, "descartes 3");
    ASSERT_TRUE(&ctx, ml_descartes_bound((long long[]){1, 0, 1}, 2) == 0, "descartes 0");

    /* S_6: dim(3,2,1) = 16; classes of S_4: (2,2)->3, (3,1)->8. */
    ASSERT_TRUE(&ctx, ml_hook_dimension((int[]){3, 2, 1}, 3) == 16, "hook 321");
    ASSERT_TRUE(&ctx, ml_hook_dimension((int[]){6}, 1) == 1, "hook trivial");
    ASSERT_TRUE(&ctx, ml_hook_dimension((int[]){1, 2}, 2) == 0, "hook non-ferrers");
    ASSERT_TRUE(&ctx, ml_sn_class_size((int[]){2, 2}, 2) == 3, "class 2+2");
    ASSERT_TRUE(&ctx, ml_sn_class_size((int[]){3, 1}, 2) == 8, "class 3+1");
    ASSERT_TRUE(&ctx, ml_sn_class_size((int[]){1, 1, 1, 1}, 4) == 1, "class identity");

    /* Rodrigues: 90 deg about z. */
    {
        double R[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double ax[] = {0, 0, 1};
        ml_status_t st = ml_rodrigues_exp(ax, 1.5707963267948966, R);
        ASSERT_TRUE(&ctx, st == ML_SUCCESS, "rodrigues status");
        ASSERT_NEAR(&ctx, R[0], 0.0, 1e-12, "R00");
        ASSERT_NEAR(&ctx, R[1], -1.0, 1e-12, "R01");
        ASSERT_NEAR(&ctx, R[3], 1.0, 1e-12, "R10");
        ASSERT_NEAR(&ctx, R[8], 1.0, 1e-12, "R22");
        ASSERT_TRUE(&ctx, ml_rodrigues_exp((double[]){0, 0, 0}, 1.0, R) == ML_ERR_INVALID_ARG, "rodrigues zero axis");
    }

    /* Ei / li / psi. */
    ASSERT_NEAR(&ctx, ml_ei(1.0), 1.8951178163559368, 1e-9, "Ei(1)");
    ASSERT_NEAR(&ctx, ml_ei(-1.0), -0.21938393439552062, 1e-9, "Ei(-1)");
    ASSERT_TRUE(&ctx, ml_isinf(ml_ei(0.0)) && ml_ei(0.0) < 0.0, "Ei(0)=-inf");
    ASSERT_NEAR(&ctx, ml_li(1000000.0), 78627.5492, 0.01, "li(1e6) PNT");
    ASSERT_TRUE(&ctx, ml_isnan(ml_li(1.0)), "li(1) NaN");
    ASSERT_NEAR(&ctx, ml_psi_chebyshev(100.0), 94.045, 0.01, "psi(100)");
    ASSERT_NEAR(&ctx, ml_psi_chebyshev(10.0), 7.832, 0.01, "psi(10)");

    return ml_test_summary(&ctx);
}
