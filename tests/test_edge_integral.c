/* v11S CLOSURE IP-20: edge integral tests */
#include "test_harness.h"
#include "ml_integral.h"

int main(void) {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Edge Integral");

    ASSERT_NEAR(&ctx, ml_factorial_float(0.0), 1.0, 1e-15, "factorial_float(0)");
    ASSERT_NEAR(&ctx, ml_factorial_float(5.0), 120.0, 1e-9, "factorial_float(5)");
    ASSERT_TRUE(&ctx, ml_isnan(ml_factorial_float(-1.0)), "factorial_float negative is NaN");

    ASSERT_NEAR(&ctx, ml_gamma_new(1.0), 1.0, 1e-12, "gamma(1)");
    ASSERT_NEAR(&ctx, ml_gamma_new(2.0), 1.0, 1e-12, "gamma(2)");
    ASSERT_TRUE(&ctx, ml_isnan(ml_gamma_new(0.0)), "gamma(0) is NaN");

    double gi = ml_gamma_new(ml_make_inf(0));
    ASSERT_TRUE(&ctx, ml_isinf(gi) && gi > 0.0, "gamma(+inf) is +inf");

    ASSERT_TRUE(&ctx, ml_isnan(ml_integral_traditional(0.0, 1.0, 2.0, 0.0, 0.0)), "integral d=0 is NaN");
    ASSERT_NEAR(&ctx, ml_integral_traditional(0.0, 1.0, 2.0, 0.0, 0.0001), 1.0 / 3.0, 1e-3, "integral x^2");

    return ml_test_summary(&ctx);
}
