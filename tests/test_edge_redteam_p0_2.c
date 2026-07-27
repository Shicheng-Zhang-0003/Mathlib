/* MATHLIB_REDP2_P0_2_TEST */
/* Red Team P2 P0-2: ml_round correctness regression tests */
#include "test_harness.h"
#include "ml_core.h"

int main(void) {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Edge Red Team P0-2");

    /*
     * Basic half-away-from-zero behavior.
     */
    ASSERT_TRUE(&ctx, ml_round(0.5) == 1.0, "round(0.5) == 1");
    ASSERT_TRUE(&ctx, ml_round(-0.5) == -1.0, "round(-0.5) == -1");

    ASSERT_TRUE(&ctx, ml_round(1.5) == 2.0, "round(1.5) == 2");
    ASSERT_TRUE(&ctx, ml_round(-1.5) == -2.0, "round(-1.5) == -2");

    ASSERT_TRUE(&ctx, ml_round(2.5) == 3.0, "round(2.5) == 3");
    ASSERT_TRUE(&ctx, ml_round(-2.5) == -3.0, "round(-2.5) == -3");

    /*
     * The old x + 0.5 implementation could misround values just below 0.5.
     *
     * 0x1.fffffffffffffp-2 is the double immediately below 0.5.
     * Correct rounding must return 0, not 1.
     */
    ASSERT_TRUE(&ctx,
        ml_round(0x1.fffffffffffffp-2) == 0.0,
        "round(nextbelow(0.5)) == 0");

    double rn = ml_round(-0x1.fffffffffffffp-2);
    ASSERT_TRUE(&ctx,
        rn == 0.0 && ml_signbit(rn) != 0,
        "round(-nextbelow(0.5)) == -0");

    /*
     * Values just above 0.5 must round to 1 / -1.
     */
    ASSERT_TRUE(&ctx,
        ml_round(0x1.0000000000001p-1) == 1.0,
        "round(nextabove(0.5)) == 1");

    ASSERT_TRUE(&ctx,
        ml_round(-0x1.0000000000001p-1) == -1.0,
        "round(-nextabove(0.5)) == -1");

    /*
     * Signed zero behavior.
     */
    ASSERT_TRUE(&ctx,
        ml_round(0.0) == 0.0 && ml_signbit(ml_round(0.0)) == 0,
        "round(+0) preserves +0");

    ASSERT_TRUE(&ctx,
        ml_round(-0.0) == 0.0 && ml_signbit(ml_round(-0.0)) != 0,
        "round(-0) preserves -0");

    /*
     * Non-finite behavior.
     */
    ASSERT_TRUE(&ctx,
        ml_isnan(ml_round(ml_make_nan())),
        "round(NaN) is NaN");

    ASSERT_TRUE(&ctx,
        ml_isinf(ml_round(ml_make_inf(0))),
        "round(+inf) is +inf");

    ASSERT_TRUE(&ctx,
        ml_isinf(ml_round(-ml_make_inf(0))),
        "round(-inf) is -inf");

    /*
     * Large integers remain unchanged.
     */
    ASSERT_TRUE(&ctx,
        ml_round(1e16) == 1e16,
        "round(large integer) identity");

    ASSERT_TRUE(&ctx,
        ml_round(-1e16) == -1e16,
        "round(-large integer) identity");

    /*
     * Boundary near 2^52.
     *
     * 0x1.fffffffffffffp+51 == 2^52 - 0.5
     * It must round to 2^52.
     */
    ASSERT_TRUE(&ctx,
        ml_round(0x1.fffffffffffffp+51) == 0x1p+52,
        "round(2^52 - 0.5) == 2^52");

    ASSERT_TRUE(&ctx,
        ml_round(-0x1.fffffffffffffp+51) == -0x1p+52,
        "round(-(2^52 - 0.5)) == -2^52");

    return ml_test_summary(&ctx);
}
