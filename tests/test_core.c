#include "test_harness.h"
#include "ml_core.h"
#include "fast_math.h"
#include "ieee754.h"

int main() {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Core & IEEE-754");
    printf("=== Core & IEEE-754 Tests ===\n");

    ASSERT_TRUE(&ctx, ml_isnan(ml_make_nan()), "NaN detection");
    ASSERT_TRUE(&ctx, ml_isinf(ml_make_inf(0)), "Inf detection");
    ASSERT_TRUE(&ctx, ml_isinf(-ml_make_inf(0)), "Negative Inf detection");
    ASSERT_NEAR(&ctx, ml_fabs(-5.5), 5.5, 1e-15, "fabs");
    ASSERT_NEAR(&ctx, ml_fast_rsqrt(4.0), 0.5, 1e-4, "fast_rsqrt");

    return ml_test_summary(&ctx);
}
