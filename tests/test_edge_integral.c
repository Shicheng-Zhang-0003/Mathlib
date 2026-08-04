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

    /* MATHLIB_V12A1_GAMMA_DD2_TEST */
    /* Exact integer factorials (shortcut path) */
    ASSERT_TRUE(&ctx, ml_gamma_new(3.0) == 2.0, "gamma(3) == 2 exact");
    ASSERT_TRUE(&ctx, ml_gamma_new(4.0) == 6.0, "gamma(4) == 6 exact");
    ASSERT_TRUE(&ctx, ml_gamma_new(5.0) == 24.0, "gamma(5) == 24 exact");
    ASSERT_TRUE(&ctx, ml_gamma_new(6.0) == 120.0, "gamma(6) == 120 exact");
    ASSERT_TRUE(&ctx, ml_gamma_new(10.0) == 362880.0, "gamma(10) == 9! exact");
    ASSERT_TRUE(&ctx, ml_gamma_new(20.0) == 121645100408832000.0, "gamma(20) == 19! exact");
    ASSERT_TRUE(&ctx, ml_gamma_new(23.0) == 1.12400072777760768e21, "gamma(23) == 22! exact");

    /* Tight comparisons against mpmath oracle values (~5.4 ULP tolerance) */
    ASSERT_NEAR(&ctx, ml_gamma_new(0.5), 1.77245385090551610e+00,
                1.77245385090551610e+00 * 1.2e-15, "gamma(0.5) == sqrt(pi) tight");
    ASSERT_NEAR(&ctx, ml_gamma_new(0.1), 9.51350769866873058e+00,
                9.51350769866873058e+00 * 1.2e-15, "gamma(0.1) tight");
    ASSERT_NEAR(&ctx, ml_gamma_new(0.01), 9.94325851191506018e+01,
                9.94325851191506018e+01 * 1.2e-15, "gamma(0.01) tight");
    ASSERT_NEAR(&ctx, ml_gamma_new(0.001), 9.99423772484595474e+02,
                9.99423772484595474e+02 * 1.2e-15, "gamma(0.001) tight");
    ASSERT_NEAR(&ctx, ml_gamma_new(50.0), 6.08281864034267522e+62,
                6.08281864034267522e+62 * 1.2e-15, "gamma(50) tight");
    ASSERT_NEAR(&ctx, ml_gamma_new(100.0), 9.33262154439441533e+155,
                9.33262154439441533e+155 * 1.2e-15, "gamma(100) tight");
    ASSERT_NEAR(&ctx, ml_gamma_new(171.0), 7.25741561530799904e+306,
                7.25741561530799904e+306 * 1.2e-15, "gamma(171) tight");

    /* Reflection formula: negative non-integer arguments */
    ASSERT_NEAR(&ctx, ml_gamma_new(-0.5), -3.54490770181103221e+00,
                3.54490770181103221e+00 * 1.2e-15, "gamma(-0.5) tight");
    ASSERT_NEAR(&ctx, ml_gamma_new(-1.5), 2.36327180120735481e+00,
                2.36327180120735481e+00 * 1.2e-15, "gamma(-1.5) tight");
    ASSERT_NEAR(&ctx, ml_gamma_new(-2.5), -9.45308720482941900e-01,
                9.45308720482941900e-01 * 1.2e-15, "gamma(-2.5) tight");
    ASSERT_NEAR(&ctx, ml_gamma_new(-3.5), 2.70088205852269114e-01,
                2.70088205852269114e-01 * 1.2e-15, "gamma(-3.5) tight");
    ASSERT_TRUE(&ctx, ml_isnan(ml_gamma_new(-1.0)), "gamma(-1) pole is NaN");
    ASSERT_TRUE(&ctx, ml_isnan(ml_gamma_new(-2.0)), "gamma(-2) pole is NaN");

    /* lgamma */
    ASSERT_TRUE(&ctx, ml_lgamma(1.0) == 0.0, "lgamma(1) == 0 exact");
    ASSERT_TRUE(&ctx, ml_lgamma(2.0) == 0.0, "lgamma(2) == 0 exact");
    ASSERT_NEAR(&ctx, ml_lgamma(5.0), ml_log(24.0), 1e-12, "lgamma(5) == log(24)");
    ASSERT_NEAR(&ctx, ml_lgamma(0.5), 5.72364942924700082e-01,
                5.72364942924700082e-01 * 1.2e-15, "lgamma(0.5) tight");
    ASSERT_NEAR(&ctx, ml_lgamma(1.5), -1.20782237635245218e-01,
                1.20782237635245218e-01 * 1.2e-15, "lgamma(1.5) tight");
    ASSERT_NEAR(&ctx, ml_lgamma(50.0), 1.44565743946344895e+02,
                1.44565743946344895e+02 * 1.2e-15, "lgamma(50) tight");
    ASSERT_NEAR(&ctx, ml_lgamma(100.0), 3.59134205369575398e+02,
                3.59134205369575398e+02 * 1.2e-15, "lgamma(100) tight");
    ASSERT_NEAR(&ctx, ml_lgamma(171.0), 7.06573062245787355e+02,
                7.06573062245787355e+02 * 1.2e-15, "lgamma(171) tight");
    ASSERT_NEAR(&ctx, ml_lgamma(-1.5), 8.60047015376480983e-01,
                8.60047015376480983e-01 * 1.2e-15, "lgamma(-1.5) tight");
    ASSERT_NEAR(&ctx, ml_lgamma(-2.5), -5.62437164976740539e-02,
                5.62437164976740539e-02 * 1.2e-15, "lgamma(-2.5) tight");
    ASSERT_TRUE(&ctx, ml_isinf(ml_lgamma(0.0)), "lgamma(0) is inf");
    ASSERT_TRUE(&ctx, ml_isnan(ml_lgamma(ml_make_nan())), "lgamma(NaN) is NaN");

    return ml_test_summary(&ctx);
}
