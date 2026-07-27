/* v11S AUDIT IP-2 regression tests */
#include "test_harness.h"
#include "ml_trig.h"
#include "ml_core.h"

int main(void) {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Edge Audit IP-2");

    /* NaN propagation */
    ASSERT_TRUE(&ctx, ml_isnan(ml_atan2(ml_make_nan(), 0.0)), "atan2(NaN,0) is NaN");
    ASSERT_TRUE(&ctx, ml_isnan(ml_atan2(0.0, ml_make_nan())), "atan2(0,NaN) is NaN");
    ASSERT_TRUE(&ctx, ml_isnan(ml_atan2(ml_make_nan(), ml_make_inf(0))), "atan2(NaN,inf) is NaN");
    ASSERT_TRUE(&ctx, ml_isnan(ml_atan2(ml_make_inf(0), ml_make_nan())), "atan2(inf,NaN) is NaN");
    ASSERT_TRUE(&ctx, ml_isnan(ml_atan2(ml_make_nan(), ml_make_nan())), "atan2(NaN,NaN) is NaN");

    /* Signed zero / zero-zero cases */
    ASSERT_TRUE(&ctx, ml_signbit(ml_atan2(-0.0, 1.0)) != 0, "atan2(-0,+1) preserves -0");
    ASSERT_TRUE(&ctx, ml_signbit(ml_atan2(0.0, 1.0)) == 0, "atan2(+0,+1) preserves +0");
    ASSERT_TRUE(&ctx, ml_signbit(ml_atan2(-0.0, ml_make_inf(0))) != 0, "atan2(-0,+inf) preserves -0");
    ASSERT_TRUE(&ctx, ml_signbit(ml_atan2(0.0, ml_make_inf(0))) == 0, "atan2(+0,+inf) preserves +0");

    ASSERT_NEAR(&ctx, ml_atan2(0.0, -1.0), ML_PI, 1e-15, "atan2(+0,-1) == +pi");
    ASSERT_NEAR(&ctx, ml_atan2(-0.0, -1.0), -ML_PI, 1e-15, "atan2(-0,-1) == -pi");

    ASSERT_NEAR(&ctx, ml_atan2(0.0, -0.0), ML_PI, 1e-15, "atan2(+0,-0) == +pi");
    ASSERT_NEAR(&ctx, ml_atan2(-0.0, -0.0), -ML_PI, 1e-15, "atan2(-0,-0) == -pi");
    ASSERT_TRUE(&ctx, ml_signbit(ml_atan2(-0.0, 0.0)) != 0, "atan2(-0,+0) preserves -0");
    ASSERT_TRUE(&ctx, ml_signbit(ml_atan2(0.0, 0.0)) == 0, "atan2(+0,+0) preserves +0");

    /* Infinity behavior */
    ASSERT_TRUE(&ctx, ml_signbit(ml_atan2(-1.0, ml_make_inf(0))) != 0, "atan2(-finite,+inf) preserves -0");
    ASSERT_TRUE(&ctx, ml_signbit(ml_atan2(1.0, ml_make_inf(0))) == 0, "atan2(+finite,+inf) preserves +0");

    ASSERT_NEAR(&ctx, ml_atan2(1.0, -ml_make_inf(0)), ML_PI, 1e-15, "atan2(+finite,-inf) == +pi");
    ASSERT_NEAR(&ctx, ml_atan2(-1.0, -ml_make_inf(0)), -ML_PI, 1e-15, "atan2(-finite,-inf) == -pi");

    ASSERT_NEAR(&ctx, ml_atan2(ml_make_inf(0), 1.0), ML_PI / 2.0, 1e-15, "atan2(+inf,finite) == +pi/2");
    ASSERT_NEAR(&ctx, ml_atan2(-ml_make_inf(0), -1.0), -ML_PI / 2.0, 1e-15, "atan2(-inf,finite) == -pi/2");

    ASSERT_NEAR(&ctx, ml_atan2(ml_make_inf(0), ml_make_inf(0)), ML_PI / 4.0, 1e-15, "atan2(+inf,+inf) == +pi/4");
    ASSERT_NEAR(&ctx, ml_atan2(ml_make_inf(0), -ml_make_inf(0)), 3.0 * ML_PI / 4.0, 1e-15, "atan2(+inf,-inf) == +3pi/4");
    ASSERT_NEAR(&ctx, ml_atan2(-ml_make_inf(0), -ml_make_inf(0)), -3.0 * ML_PI / 4.0, 1e-15, "atan2(-inf,-inf) == -3pi/4");
    ASSERT_NEAR(&ctx, ml_atan2(-ml_make_inf(0), ml_make_inf(0)), -ML_PI / 4.0, 1e-15, "atan2(-inf,+inf) == -pi/4");

    ASSERT_NEAR(&ctx, ml_atan2(ml_make_inf(0), 0.0), ML_PI / 2.0, 1e-15, "atan2(+inf,+0) == +pi/2");
    ASSERT_NEAR(&ctx, ml_atan2(-ml_make_inf(0), -0.0), -ML_PI / 2.0, 1e-15, "atan2(-inf,-0) == -pi/2");

    ASSERT_NEAR(&ctx, ml_atan2(0.0, -ml_make_inf(0)), ML_PI, 1e-15, "atan2(+0,-inf) == +pi");
    ASSERT_NEAR(&ctx, ml_atan2(-0.0, -ml_make_inf(0)), -ML_PI, 1e-15, "atan2(-0,-inf) == -pi");

    return ml_test_summary(&ctx);
}
