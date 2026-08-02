/* v11S CLOSURE IP-20: edge integral tests */
#include "test_harness.h"
#include "ml_integral.h"

int main(void) {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Edge Integral");

    ASSERT_NEAR(&ctx, ml_factorial_float(0.0), 1.0, 1e-15, "factorial_float(0)");
    ASSERT_NEAR(&ctx, ml_factorial_float(5.0), 120.0, 1e-9, "factorial_float(5)");
    ASSERT_TRUE(&ctx, ml_isnan(ml_factorial_float(-1.0)), "factorial_float negative is NaN");

    ASSERT_NEAR(&ctx, ml_gamma_new(1.0), 1.0, 1e-14, "gamma(1)");
    ASSERT_NEAR(&ctx, ml_gamma_new(2.0), 1.0, 1e-14, "gamma(2)");
    ASSERT_TRUE(&ctx, ml_isnan(ml_gamma_new(0.0)), "gamma(0) is NaN");

    double gi = ml_gamma_new(ml_make_inf(0));
    ASSERT_TRUE(&ctx, ml_isinf(gi) && gi > 0.0, "gamma(+inf) is +inf");

    ASSERT_TRUE(&ctx, ml_isnan(ml_integral_traditional(0.0, 1.0, 2.0, 0.0, 0.0)), "integral d=0 is NaN");
    ASSERT_NEAR(&ctx, ml_integral_traditional(0.0, 1.0, 2.0, 0.0, 0.0001), 1.0 / 3.0, 1e-3, "integral x^2");


    /* MATHLIB_V12A1_GAMMA_LANCZOS_TEST */
    /* Lanczos accuracy: known exact values */
    ASSERT_NEAR(&ctx, ml_gamma_new(3.0), 2.0, 1e-13, "gamma(3) == 2");
    ASSERT_NEAR(&ctx, ml_gamma_new(4.0), 6.0, 1e-12, "gamma(4) == 6");
    ASSERT_NEAR(&ctx, ml_gamma_new(5.0), 24.0, 1e-11, "gamma(5) == 24");
    ASSERT_NEAR(&ctx, ml_gamma_new(6.0), 120.0, 1e-10, "gamma(6) == 120");
    ASSERT_NEAR(&ctx, ml_gamma_new(0.5), 1.7724538509055159, 1e-13, "gamma(0.5) == sqrt(pi)");

    /* Reflection formula: negative non-integer arguments */
    ASSERT_NEAR(&ctx, ml_gamma_new(-0.5), -3.5449077018110318, 1e-12, "gamma(-0.5)");
    ASSERT_NEAR(&ctx, ml_gamma_new(-1.5), 2.3632718012073548, 1e-12, "gamma(-1.5)");
    ASSERT_TRUE(&ctx, ml_isnan(ml_gamma_new(-1.0)), "gamma(-1) pole is NaN");
    ASSERT_TRUE(&ctx, ml_isnan(ml_gamma_new(-2.0)), "gamma(-2) pole is NaN");

    /* lgamma basic checks */
    ASSERT_NEAR(&ctx, ml_lgamma(1.0), 0.0, 1e-14, "lgamma(1) == 0");
    ASSERT_NEAR(&ctx, ml_lgamma(2.0), 0.0, 1e-14, "lgamma(2) == 0");
    ASSERT_NEAR(&ctx, ml_lgamma(5.0), ml_log(24.0), 1e-12, "lgamma(5) == log(24)");
    ASSERT_TRUE(&ctx, ml_isinf(ml_lgamma(0.0)), "lgamma(0) is inf");
    ASSERT_TRUE(&ctx, ml_isnan(ml_lgamma(ml_make_nan())), "lgamma(NaN) is NaN");

    return ml_test_summary(&ctx);
}
