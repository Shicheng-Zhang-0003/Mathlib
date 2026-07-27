/* MATHLIB_CLOSURE_P2_P0_5_TEST */
/* Red Team P2 P0-5: CORDIC non-finite regression tests */
#include "test_harness.h"
#include "ml_core.h"
#include "internal/cordic.h"

int main(void) {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Edge Red Team P0-5");

    double s, c;

    ml_cordic_sincos(ml_make_nan(), &s, &c);
    ASSERT_TRUE(&ctx,
        ml_isnan(s) && ml_isnan(c),
        "CORDIC NaN propagates");

    ml_cordic_sincos(ml_make_inf(0), &s, &c);
    ASSERT_TRUE(&ctx,
        ml_isnan(s) && ml_isnan(c),
        "CORDIC +inf propagates");

    ml_cordic_sincos(-ml_make_inf(0), &s, &c);
    ASSERT_TRUE(&ctx,
        ml_isnan(s) && ml_isnan(c),
        "CORDIC -inf propagates");

    ml_cordic_sincos(0.0, &s, &c);
    ASSERT_NEAR(&ctx, s, 0.0, 1e-6, "CORDIC sin(0)");
    ASSERT_NEAR(&ctx, c, 1.0, 1e-6, "CORDIC cos(0)");

    ml_cordic_sincos(ML_PI / 2.0, &s, &c);
    ASSERT_NEAR(&ctx, s, 1.0, 1e-5, "CORDIC sin(pi/2)");
    ASSERT_NEAR(&ctx, ml_fabs(c), 0.0, 1e-5, "CORDIC cos(pi/2) near zero");

    return ml_test_summary(&ctx);
}
