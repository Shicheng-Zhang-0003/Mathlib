/* v11S CLOSURE IP-20: edge trig tests */
#include "test_harness.h"
#include "ml_trig.h"

int main(void) {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Edge Trig");

    ASSERT_TRUE(&ctx, ml_isnan(ml_sin(ml_make_nan())), "sin(NaN) is NaN");
    ASSERT_TRUE(&ctx, ml_isnan(ml_sin(ml_make_inf(0))), "sin(Inf) is NaN");
    ASSERT_TRUE(&ctx, ml_isnan(ml_cos(ml_make_nan())), "cos(NaN) is NaN");
    ASSERT_TRUE(&ctx, ml_isnan(ml_cos(ml_make_inf(0))), "cos(Inf) is NaN");

    ASSERT_NEAR(&ctx, ml_tan(0.0), 0.0, 1e-15, "tan(0)");
    ASSERT_NEAR(&ctx, ml_tan(ML_PI / 4.0), 1.0, 1e-14, "tan(pi/4)");

    ASSERT_NEAR(&ctx, ml_atan(ml_make_inf(0)), ML_PI / 2.0, 1e-15, "atan(+inf)");
    ASSERT_NEAR(&ctx, ml_atan(-ml_make_inf(0)), -ML_PI / 2.0, 1e-15, "atan(-inf)");
    ASSERT_TRUE(&ctx, ml_isnan(ml_atan(ml_make_nan())), "atan(NaN) is NaN");

    ASSERT_TRUE(&ctx, ml_isnan(ml_asin(2.0)), "asin domain");
    ASSERT_TRUE(&ctx, ml_isnan(ml_acos(2.0)), "acos domain");

    ASSERT_NEAR(&ctx, ml_atan2(0.0, -0.0), ML_PI, 1e-14, "atan2(+0,-0)");
    ASSERT_NEAR(&ctx, ml_atan2(-0.0, -0.0), -ML_PI, 1e-14, "atan2(-0,-0)");
    ASSERT_NEAR(&ctx, ml_atan2(0.0, 0.0), 0.0, 1e-14, "atan2(+0,+0)");

    ASSERT_TRUE(&ctx, !ml_isnan(ml_sin(1e10)) && !ml_isinf(ml_sin(1e10)), "sin(1e10) finite");
    ASSERT_TRUE(&ctx, ml_isnan(ml_sin(1e50)), "sin(1e50) safely NaN");

    return ml_test_summary(&ctx);
}
