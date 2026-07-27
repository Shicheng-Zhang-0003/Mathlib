/* v11S CLOSURE IP-20: edge numerical tests */
#include "test_harness.h"
#include "ml_numerical.h"

static double f_quad(double x) {
    return x * x - 4.0;
}

static double df_quad(double x) {
    return 2.0 * x;
}

int main(void) {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Edge Numerical");

    ASSERT_TRUE(&ctx, ml_isnan(ml_newton_raphson(NULL, df_quad, 3.0, 1e-12, 100)), "newton rejects NULL f");
    ASSERT_TRUE(&ctx, ml_isnan(ml_newton_raphson(f_quad, NULL, 3.0, 1e-12, 100)), "newton rejects NULL df");

    ASSERT_TRUE(&ctx, ml_isnan(ml_bisection(NULL, 1.0, 3.0, 1e-12, 100)), "bisection rejects NULL f");
    ASSERT_TRUE(&ctx, ml_isnan(ml_bisection(f_quad, 3.0, 4.0, 1e-12, 100)), "bisection rejects same-sign endpoints");

    ASSERT_TRUE(&ctx, ml_isnan(ml_derivative(NULL, 1.0, 1e-5)), "derivative rejects NULL f");
    ASSERT_TRUE(&ctx, ml_isnan(ml_derivative(f_quad, 1.0, 0.0)), "derivative rejects h=0");

    ASSERT_TRUE(&ctx, ml_isnan(ml_integral_simpson(NULL, 0.0, 1.0, 2)), "simpson rejects NULL f");
    ASSERT_TRUE(&ctx, ml_isnan(ml_integral_simpson(f_quad, 0.0, 1.0, 3)), "simpson rejects odd n");

    ASSERT_NEAR(&ctx, ml_integral_simpson(f_quad, 0.0, 2.0, 100), -16.0 / 3.0, 1e-6, "simpson integral value");

    return ml_test_summary(&ctx);
}
