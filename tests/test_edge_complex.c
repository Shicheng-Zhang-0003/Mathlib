/* v11S CLOSURE IP-20: edge complex tests */
#include "test_harness.h"
#include "ml_complex.h"

int main(void) {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Edge Complex");

    cplx big = {1e308, 1e308};
    double abs_big = ml_cplx_abs(big);

    ASSERT_TRUE(&ctx, ml_isfinite(abs_big), "complex abs overflow-safe");
    ASSERT_NEAR(&ctx, abs_big, 1.4142135623730951e308, 1e296, "complex abs value");

    cplx zero = {0.0, 0.0};
    ASSERT_NEAR(&ctx, ml_cplx_abs(zero), 0.0, 0.0, "complex abs zero");

    ASSERT_NEAR(&ctx, ml_cplx_arg((cplx){1.0, 0.0}), 0.0, 1e-15, "arg(1)");
    ASSERT_NEAR(&ctx, ml_cplx_arg((cplx){0.0, 1.0}), ML_PI / 2.0, 1e-15, "arg(i)");
    ASSERT_NEAR(&ctx, ml_cplx_arg((cplx){-1.0, 0.0}), ML_PI, 1e-15, "arg(-1)");
    ASSERT_NEAR(&ctx, ml_cplx_arg((cplx){0.0, -1.0}), -ML_PI / 2.0, 1e-15, "arg(-i)");

    ASSERT_TRUE(&ctx, ml_isnan(ml_cplx_arg((cplx){ml_make_nan(), 0.0})), "arg(NaN) is NaN");

    cplx div_zero = ml_cplx_div((cplx){1.0, 0.0}, (cplx){0.0, 0.0});
    ASSERT_TRUE(&ctx, ml_isnan(div_zero.real) && ml_isnan(div_zero.imag), "complex div by zero gives NaN");

    cplx exp_nan = ml_cplx_exponential((cplx){ml_make_nan(), 0.0});
    ASSERT_TRUE(&ctx, ml_isnan(exp_nan.real), "complex exp NaN guard");

    cplx log_zero = ml_cplx_logarithm((cplx){0.0, 0.0});
    ASSERT_TRUE(&ctx, ml_isinf(log_zero.real) && log_zero.real < 0.0, "complex log(0) real is -inf");

    return ml_test_summary(&ctx);
}
