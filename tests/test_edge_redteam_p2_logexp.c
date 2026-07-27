/* MATHLIB_CLOSURE_P2_LOG_EXP_ROUNDTRIP_TEST */
/* Red Team P2: log(exp(x)) roundtrip regression tests */
#include "test_harness.h"
#include "ml_core.h"
#include "ml_exp_log.h"

#include <stdint.h>
#include <string.h>

static uint64_t redteam_ulp_distance(double a, double b) {
    uint64_t ia, ib;

    memcpy(&ia, &a, sizeof(uint64_t));
    memcpy(&ib, &b, sizeof(uint64_t));

    if (ia >> 63) ia = 0x8000000000000000ULL - ia;
    if (ib >> 63) ib = 0x8000000000000000ULL - ib;

    return ia > ib ? ia - ib : ib - ia;
}

static void check_roundtrip(ml_test_ctx_t *ctx, double x) {
    double e = ml_exp(x);
    double l = ml_log(e);

    ASSERT_TRUE(ctx,
        ml_isfinite(e) && ml_isfinite(l),
        "exp/log roundtrip finite");

    if (!ml_isfinite(e) || !ml_isfinite(l)) {
        return;
    }

    double diff = ml_fabs(l - x);
    uint64_t ulps = redteam_ulp_distance(l, x);

    ASSERT_TRUE(ctx,
        diff <= 1e-13 || ulps <= 8,
        "log(exp(x)) roundtrip within 1e-13 absolute or 8 ULP");
}

int main(void) {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Edge Red Team P2 Log/Exp Roundtrip");

    /*
     * These values came from deterministic fuzz_god_mode seed 123456789.
     * They are in the double range where 1 ULP is about 1.14e-13.
     */
    check_roundtrip(&ctx, 586.978120536999);
    check_roundtrip(&ctx, 627.643520770428);
    check_roundtrip(&ctx, 671.779855933938);
    check_roundtrip(&ctx, 674.258560721883);
    check_roundtrip(&ctx, 628.762675276381);
    check_roundtrip(&ctx, 588.572237449033);
    check_roundtrip(&ctx, 609.755943347586);
    check_roundtrip(&ctx, 605.736818912316);
    check_roundtrip(&ctx, 688.953238394555);
    check_roundtrip(&ctx, 628.651212727954);
    check_roundtrip(&ctx, 651.39389487514);
    check_roundtrip(&ctx, 672.896168042392);
    check_roundtrip(&ctx, 630.227606571385);
    check_roundtrip(&ctx, 632.544160650434);
    check_roundtrip(&ctx, 654.068812567841);
    check_roundtrip(&ctx, 628.888303241175);
    check_roundtrip(&ctx, 544.08798345555);

    return ml_test_summary(&ctx);
}
