/* MATHLIB_REDP2_P0_1_TEST */
/* Red Team P2 P0-1: ldexp subnormal overflow regression tests */
#include "test_harness.h"
#include "ml_core.h"

#include <stdint.h>
#include <string.h>

int main(void) {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Edge Red Team P0-1");

    uint64_t min_bits = 1;
    double min_sub;
    memcpy(&min_sub, &min_bits, sizeof(double));

    /*
     * Subnormal + large positive exponent must not falsely overflow.
     */
    double r2046 = ml_ldexp_pure(min_sub, 2046);
    ASSERT_TRUE(&ctx,
        ml_isfinite(r2046),
        "ldexp(min_subnormal, 2046) is finite");

    ASSERT_TRUE(&ctx,
        r2046 == 0x1p972,
        "ldexp(min_subnormal, 2046) == 2^972");

    double rneg = ml_ldexp_pure(-min_sub, 2046);
    ASSERT_TRUE(&ctx,
        ml_isfinite(rneg) && ml_signbit(rneg) != 0,
        "ldexp(-min_subnormal, 2046) preserves sign");

    ASSERT_TRUE(&ctx,
        rneg == -0x1p972,
        "ldexp(-min_subnormal, 2046) == -2^972");

    /*
     * Finite boundary near the top of the double range.
     */
    double r2097 = ml_ldexp_pure(min_sub, 2097);
    ASSERT_TRUE(&ctx,
        ml_isfinite(r2097),
        "ldexp(min_subnormal, 2097) is finite");

    ASSERT_TRUE(&ctx,
        r2097 == 0x1p1023,
        "ldexp(min_subnormal, 2097) == 2^1023");

    double r2098 = ml_ldexp_pure(min_sub, 2098);
    ASSERT_TRUE(&ctx,
        ml_isinf(r2098) && r2098 > 0.0,
        "ldexp(min_subnormal, 2098) overflows to +inf");

    /*
     * Normal-input sanity boundaries.
     */
    double n1023 = ml_ldexp_pure(1.0, 1023);
    ASSERT_TRUE(&ctx,
        ml_isfinite(n1023),
        "ldexp(1, 1023) is finite");

    ASSERT_TRUE(&ctx,
        n1023 == 0x1p1023,
        "ldexp(1, 1023) == 2^1023");

    double p1024 = ml_ldexp_pure(1.0, 1024);
    ASSERT_TRUE(&ctx,
        ml_isinf(p1024) && p1024 > 0.0,
        "ldexp(1, 1024) overflows to +inf");

    /*
     * Existing underflow behavior must remain intact.
     */
    ASSERT_TRUE(&ctx,
        ml_ldexp_pure(1.0, -1074) == min_sub,
        "ldexp(1, -1074) == min_subnormal");

    ASSERT_TRUE(&ctx,
        ml_ldexp_pure(1.0, -1075) == 0.0,
        "ldexp(1, -1075) underflows to zero");

    ASSERT_TRUE(&ctx,
        ml_signbit(ml_ldexp_pure(-1.0, -1075)) != 0,
        "ldexp(-1, -1075) preserves negative zero sign");

    return ml_test_summary(&ctx);
}
