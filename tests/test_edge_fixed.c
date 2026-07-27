/* v11S CLOSURE IP-20: edge fixed-point tests */
#include "test_harness.h"
#include "ml_fixed_point.h"

int main(void) {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Edge Fixed-Point");

    ASSERT_TRUE(&ctx, ml_fixed_mul(2147483647, 2147483647) == 2147483647, "fixed mul saturates max");
    ASSERT_TRUE(&ctx, ml_fixed_mul((-2147483647 - 1), 2147483647) == (-2147483647 - 1), "fixed mul saturates min");
    ASSERT_TRUE(&ctx, ml_fixed_mul(2 << 16, 3 << 16) == (6 << 16), "fixed mul 2*3");

    ml_q16_16_t s = 0;
    ml_q16_16_t c = 0;

    ml_cordic_sincos_fixed(0, &s, &c);
    ASSERT_NEAR(&ctx, (double)s / 65536.0, 0.0, 1e-3, "fixed sin(0)");
    ASSERT_NEAR(&ctx, (double)c / 65536.0, 1.0, 1e-3, "fixed cos(0)");

    ml_cordic_sincos_fixed(ML_FIXED_HALF_PI, &s, &c);
    ASSERT_NEAR(&ctx, (double)s / 65536.0, 1.0, 1e-2, "fixed sin(pi/2)");
    ASSERT_NEAR(&ctx, (double)c / 65536.0, 0.0, 1e-2, "fixed cos(pi/2)");

    return ml_test_summary(&ctx);
}
