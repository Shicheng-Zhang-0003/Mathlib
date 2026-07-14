#include "test_harness.h"
#include "ml_trig.h"
#include "ml_exp_log.h"

int main() {
    printf("=== Trig & Exp/Log Tests ===\n");
    double pi = 3.14159265358979323846;
    ASSERT_NEAR(ml_sin(pi/2.0), 1.0, 1e-14, "sin(pi/2)");
    ASSERT_NEAR(ml_cos(pi), -1.0, 1e-14, "cos(pi)");
    ASSERT_NEAR(ml_log(1.0), 0.0, 1e-15, "log(1)");
    ASSERT_NEAR(ml_exp(0.0), 1.0, 1e-15, "exp(0)");
    return test_harness_summary("Trig/Exp");
}