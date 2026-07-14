#include "test_harness.h"
#include "ml_core.h"
#include "fast_math.h"
#include "ieee754.h"

int main() {
    printf("=== Core & IEEE-754 Tests ===\n");
    ASSERT_TRUE(ml_isnan(0.0/0.0), "NaN detection");
    ASSERT_TRUE(ml_isinf(1.0/0.0), "Inf detection");
    ASSERT_NEAR(ml_fabs(-5.5), 5.5, 1e-15, "fabs");
    ASSERT_NEAR(ml_fast_rsqrt(4.0), 0.5, 1e-4, "fast_rsqrt");
    return test_harness_summary("Core");
}