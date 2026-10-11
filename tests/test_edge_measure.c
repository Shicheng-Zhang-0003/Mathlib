/* Doctorate measure pins: Dirichlet, simple sums, Lp, Holder/Minkowski. */
#include "test_harness.h"
#include "ml_measure.h"
#include "ml_core.h"

int main(void) {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Edge Measure");

    ASSERT_NEAR(&ctx, ml_lebesgue_dirichlet(0, 1), 0.0, 0.0, "dirichlet [0,1]");
    ASSERT_NEAR(&ctx, ml_lebesgue_dirichlet(2, 2), 0.0, 0.0, "dirichlet point");
    ASSERT_TRUE(&ctx, ml_isnan(ml_lebesgue_dirichlet(1, 0)), "dirichlet a>b NaN");
    ASSERT_TRUE(&ctx, ml_isnan(ml_lebesgue_dirichlet(0, ml_make_inf(0))), "dirichlet inf NaN");
    {
        double v[] = {2.0, 3.0}, m[] = {0.5, 0.25};
        ASSERT_NEAR(&ctx, ml_simple_integral(v, m, 2), 1.75, 1e-15, "simple integral");
        ASSERT_TRUE(&ctx, ml_isnan(ml_simple_integral(v, (double[]){-1, 1}, 2)), "negative measure NaN");
    }
    {
        /* ||1||_2 on 4 unit samples, dx=0.25 -> 1. */
        double f[] = {1, 1, 1, 1};
        ASSERT_NEAR(&ctx, ml_lp_norm_data(f, 4, 0.25, 2.0), 1.0, 1e-12, "lp const");
        ASSERT_NEAR(&ctx, ml_seq_lp_norm((double[]){3, 4}, 2, 2.0), 5.0, 1e-12, "l2 3-4-5");
        ASSERT_TRUE(&ctx, ml_isnan(ml_lp_norm_data(f, 4, 0.25, 0.5)), "p<1 NaN");
        ASSERT_TRUE(&ctx, ml_isnan(ml_lp_norm_data(NULL, 4, 0.25, 2.0)), "null NaN");
    }
    {
        double f[] = {1, 2, 3, 4}, g[] = {4, 3, 2, 1};
        ASSERT_TRUE(&ctx, ml_holder_check(f, g, 4, 0.25, 2.0, 1e-9) == 1, "holder holds");
        ASSERT_TRUE(&ctx, ml_minkowski_check(f, g, 4, 0.25, 2.0, 1e-9) == 1, "minkowski holds");
        ASSERT_TRUE(&ctx, ml_holder_check(f, NULL, 4, 0.25, 2.0, 1e-9) == 0, "holder null");
        ASSERT_TRUE(&ctx, ml_minkowski_check(f, g, 4, 0.25, 1.0, 1e-9) == 1, "minkowski p=1");
    }
    ASSERT_NEAR(&ctx, ml_markov_bound(2.0, 4.0), 0.5, 1e-15, "markov bound");
    ASSERT_TRUE(&ctx, ml_isnan(ml_markov_bound(2.0, 0.0)), "markov t=0 NaN");

    return ml_test_summary(&ctx);
}
