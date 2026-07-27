/* v11S CLOSURE IP-20: edge core tests */
#include "test_harness.h"
#include "ml_core.h"
#include <math.h>
#include <stdint.h>
#include <string.h>

int main(void) {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Edge Core");

    uint64_t min_bits = 1;
    double min_subnormal;
    memcpy(&min_subnormal, &min_bits, sizeof(double));

    ASSERT_TRUE(&ctx, ml_signbit(ml_sqrt(-0.0)) != 0, "sqrt(-0) preserves sign");
    ASSERT_TRUE(&ctx, ml_signbit(ml_round(-0.0)) != 0, "round(-0) preserves sign");
    ASSERT_TRUE(&ctx, ml_signbit(ml_fabs(-0.0)) == 0, "fabs(-0) gives +0");

    ASSERT_TRUE(&ctx, ml_ldexp_pure(1.0, -1074) == min_subnormal, "ldexp min subnormal");
    ASSERT_TRUE(&ctx, ml_ldexp_pure(1.0, -1075) == 0.0, "ldexp underflow to zero");
    ASSERT_TRUE(&ctx, ml_signbit(ml_ldexp_pure(-1.0, -1075)) != 0, "ldexp underflow preserves sign");

    int e = 0;
    double m = ml_frexp_pure(min_subnormal, &e);
    ASSERT_NEAR(&ctx, m, 0.5, 0.0, "frexp min subnormal mantissa");
    ASSERT_TRUE(&ctx, e == -1073, "frexp min subnormal exponent");

    ASSERT_NEAR(&ctx, ml_fmod(5.0, 3.0), 2.0, 1e-15, "fmod(5,3)");
    ASSERT_NEAR(&ctx, ml_fmod(-5.0, 3.0), -2.0, 1e-15, "fmod(-5,3)");
    ASSERT_NEAR(&ctx, ml_fmod(5.0, -3.0), 2.0, 1e-15, "fmod(5,-3)");

    double x = 7.0 * min_subnormal;
    double y = 2.0 * min_subnormal;
    ASSERT_TRUE(&ctx, ml_fmod(x, y) == min_subnormal, "fmod subnormal exact");

    ASSERT_TRUE(&ctx, ml_signbit(ml_fmod(-0.0, 1.0)) != 0, "fmod(-0,1) preserves sign");
    ASSERT_TRUE(&ctx, ml_isnan(ml_fmod(1.0, 0.0)), "fmod(x,0) is NaN");
    ASSERT_TRUE(&ctx, ml_isnan(ml_fmod(ml_make_inf(0), 1.0)), "fmod(inf,finite) is NaN");
    ASSERT_TRUE(&ctx, ml_fmod(1.0, ml_make_inf(0)) == 1.0, "fmod(finite,inf) is x");

    double bigx = 1.23456789e+300;
    double bigy = 9.87654321e-300;
    double got = ml_fmod(bigx, bigy);
    double expected = fmod(bigx, bigy);
    ASSERT_NEAR(&ctx, got, expected, ml_fabs(expected) * 1e-13 + 1e-315, "fmod large quotient vs libc");

    return ml_test_summary(&ctx);
}
