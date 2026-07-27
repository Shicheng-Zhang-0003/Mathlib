/* v11S CLOSURE IP-20: edge hyperbolic tests */
#include "test_harness.h"
#include "ml_exp_log.h"

int main(void) {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Edge Hyperbolic");

    ASSERT_TRUE(&ctx, ml_tanh(1000.0) == 1.0, "tanh(+large) == 1");
    ASSERT_TRUE(&ctx, ml_tanh(-1000.0) == -1.0, "tanh(-large) == -1");
    ASSERT_TRUE(&ctx, ml_signbit(ml_tanh(-0.0)) != 0, "tanh(-0) preserves sign");

    ASSERT_TRUE(&ctx, ml_isinf(ml_sinh(1000.0)) && ml_sinh(1000.0) > 0.0, "sinh(+large) == +inf");
    ASSERT_TRUE(&ctx, ml_isinf(ml_sinh(-1000.0)) && ml_sinh(-1000.0) < 0.0, "sinh(-large) == -inf");
    ASSERT_TRUE(&ctx, ml_isinf(ml_cosh(1000.0)) && ml_cosh(1000.0) > 0.0, "cosh(large) == +inf");

    ASSERT_TRUE(&ctx, ml_isfinite(ml_asinh(1e308)), "asinh(large) finite");
    ASSERT_TRUE(&ctx, ml_isfinite(ml_acosh(1e308)), "acosh(large) finite");

    ASSERT_TRUE(&ctx, ml_isnan(ml_acosh(0.5)), "acosh(<1) is NaN");
    ASSERT_TRUE(&ctx, ml_isnan(ml_atanh(1.0)), "atanh(1) is NaN");
    ASSERT_TRUE(&ctx, ml_isnan(ml_atanh(-1.0)), "atanh(-1) is NaN");
    ASSERT_NEAR(&ctx, ml_atanh(0.0), 0.0, 0.0, "atanh(0) == 0");

    return ml_test_summary(&ctx);
}
