#include "test_harness.h"
#include "trigonometry.h"
#include "exponential.h"

int main() {
    printf("=== Trig & Exp/Log Tests ===\n");
    double pi = 3.14159265358979323846;
    ASSERT_NEAR(sine(pi/2.0), 1.0, 1e-14, "sin(pi/2)");
    ASSERT_NEAR(cosine(pi), -1.0, 1e-14, "cos(pi)");
    ASSERT_NEAR(logarithm(1.0), 0.0, 1e-15, "log(1)");
    ASSERT_NEAR(exponential(0.0), 1.0, 1e-15, "exp(0)");
    return test_harness_summary("Trig/Exp");
}