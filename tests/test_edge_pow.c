/* v11S CLOSURE IP-20: edge pow tests */
#include "test_harness.h"
#include "ml_exp_log.h"

int main(void) {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Edge Pow");

    ASSERT_NEAR(&ctx, ml_pow(2.0, 3.0), 8.0, 1e-14, "pow(2,3)");
    ASSERT_NEAR(&ctx, ml_pow(-2.0, 3.0), -8.0, 1e-14, "pow(-2,3)");
    ASSERT_NEAR(&ctx, ml_pow(-2.0, 2.0), 4.0, 1e-14, "pow(-2,2)");

    ASSERT_TRUE(&ctx, ml_pow(0.0, 0.0) == 1.0, "pow(0,0) == 1");
    ASSERT_TRUE(&ctx, ml_isinf(ml_pow(0.0, -1.0)) && ml_pow(0.0, -1.0) > 0.0, "pow(+0,-1) == +inf");
    ASSERT_TRUE(&ctx, ml_isinf(ml_pow(-0.0, -1.0)) && ml_signbit(ml_pow(-0.0, -1.0)) != 0, "pow(-0,-1) == -inf");

    ASSERT_TRUE(&ctx, ml_isnan(ml_pow(-1.0, 0.5)), "pow(negative, fractional) is NaN");
    ASSERT_TRUE(&ctx, ml_pow(1.0, ml_make_nan()) == 1.0, "pow(1,NaN) == 1");
    ASSERT_TRUE(&ctx, ml_pow(ml_make_nan(), 0.0) == 1.0, "pow(NaN,0) == 1");

    ASSERT_TRUE(&ctx, ml_isinf(ml_pow(2.0, ml_make_inf(0))) && ml_pow(2.0, ml_make_inf(0)) > 0.0, "pow(2,+inf) == +inf");
    ASSERT_TRUE(&ctx, ml_pow(2.0, -ml_make_inf(0)) == 0.0, "pow(2,-inf) == 0");
    ASSERT_TRUE(&ctx, ml_pow(0.5, ml_make_inf(0)) == 0.0, "pow(0.5,+inf) == 0");
    ASSERT_TRUE(&ctx, ml_isinf(ml_pow(0.5, -ml_make_inf(0))) && ml_pow(0.5, -ml_make_inf(0)) > 0.0, "pow(0.5,-inf) == +inf");

    return ml_test_summary(&ctx);
}
