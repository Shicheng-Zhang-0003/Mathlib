/* v11S CLOSURE IP-20: edge combinatorics tests */
#include "test_harness.h"
#include "ml_combinatorics.h"
#include <stdint.h>

int main(void) {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Edge Combinatorics");

    ASSERT_TRUE(&ctx, ml_factorial(0) == 1, "0!");
    ASSERT_TRUE(&ctx, ml_factorial(5) == 120, "5!");
    ASSERT_TRUE(&ctx, ml_factorial(20) == 2432902008176640000ULL, "20!");
    ASSERT_TRUE(&ctx, ml_factorial(21) == UINT64_MAX, "21! overflow sentinel");
    ASSERT_TRUE(&ctx, ml_factorial(-1) == UINT64_MAX, "negative factorial sentinel");

    ASSERT_TRUE(&ctx, ml_ncr(5, 2) == 10, "ncr(5,2)");
    ASSERT_TRUE(&ctx, ml_ncr(60, 30) == 118264581564861424ULL, "ncr(60,30)");
    ASSERT_TRUE(&ctx, ml_ncr(62, 31) == 465428353255261088ULL, "ncr(62,31)");
    ASSERT_TRUE(&ctx, ml_ncr(100, 50) == UINT64_MAX, "ncr overflow sentinel");
    ASSERT_TRUE(&ctx, ml_ncr(-1, 0) == 0, "ncr negative input");

    ASSERT_TRUE(&ctx, ml_npr(5, 2) == 20, "npr(5,2)");
    ASSERT_TRUE(&ctx, ml_npr(21, 21) == UINT64_MAX, "npr overflow sentinel");
    ASSERT_TRUE(&ctx, ml_npr(-1, 0) == 0, "npr negative input");

    return ml_test_summary(&ctx);
}
