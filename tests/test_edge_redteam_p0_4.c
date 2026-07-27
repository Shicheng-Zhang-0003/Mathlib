/* MATHLIB_CLOSURE_P2_P0_4_TEST */
/* Red Team P2 P0-4: hyperbolic overflow threshold regression tests */
#include "test_harness.h"
#include "ml_core.h"
#include "ml_exp_log.h"

int main(void) {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Edge Red Team P0-4");

    /*
     * Basic special-case behavior.
     */
    ASSERT_TRUE(&ctx,
        ml_isnan(ml_sinh(ml_make_nan())),
        "sinh(NaN) is NaN");

    ASSERT_TRUE(&ctx,
        ml_isnan(ml_cosh(ml_make_nan())),
        "cosh(NaN) is NaN");

    ASSERT_TRUE(&ctx,
        ml_isinf(ml_sinh(ml_make_inf(0))) && ml_sinh(ml_make_inf(0)) > 0.0,
        "sinh(+inf) is +inf");

    ASSERT_TRUE(&ctx,
        ml_isinf(ml_sinh(-ml_make_inf(0))) && ml_sinh(-ml_make_inf(0)) < 0.0,
        "sinh(-inf) is -inf");

    ASSERT_TRUE(&ctx,
        ml_isinf(ml_cosh(ml_make_inf(0))) && ml_cosh(ml_make_inf(0)) > 0.0,
        "cosh(+inf) is +inf");

    ASSERT_TRUE(&ctx,
        ml_isinf(ml_cosh(-ml_make_inf(0))) && ml_cosh(-ml_make_inf(0)) > 0.0,
        "cosh(-inf) is +inf");

    /*
     * The old threshold was approximately log(DBL_MAX):
     *
     *   709.782712893384
     *
     * But sinh/cosh overflow near:
     *
     *   log(DBL_MAX) + log(2)
     *   ~= 710.475860073944
     *
     * Therefore values around 710.0 must remain finite.
     */
    double s710 = ml_sinh(710.0);
    ASSERT_TRUE(&ctx,
        ml_isfinite(s710) && s710 > 0.0,
        "sinh(710) must be finite");

    double c710 = ml_cosh(710.0);
    ASSERT_TRUE(&ctx,
        ml_isfinite(c710) && c710 > 0.0,
        "cosh(710) must be finite");

    double s7104 = ml_sinh(710.4);
    ASSERT_TRUE(&ctx,
        ml_isfinite(s7104) && s7104 > 0.0,
        "sinh(710.4) must be finite");

    double c7104 = ml_cosh(710.4);
    ASSERT_TRUE(&ctx,
        ml_isfinite(c7104) && c7104 > 0.0,
        "cosh(710.4) must be finite");

    /*
     * Beyond the corrected hyperbolic overflow boundary, result must be inf.
     */
    ASSERT_TRUE(&ctx,
        ml_isinf(ml_sinh(710.5)) && ml_sinh(710.5) > 0.0,
        "sinh(710.5) must overflow to +inf");

    ASSERT_TRUE(&ctx,
        ml_isinf(ml_cosh(710.5)) && ml_cosh(710.5) > 0.0,
        "cosh(710.5) must overflow to +inf");

    ASSERT_TRUE(&ctx,
        ml_isinf(ml_sinh(-710.5)) && ml_sinh(-710.5) < 0.0,
        "sinh(-710.5) must overflow to -inf");

    /*
     * Existing large-argument behavior must remain intact.
     */
    ASSERT_TRUE(&ctx,
        ml_isinf(ml_sinh(1000.0)) && ml_sinh(1000.0) > 0.0,
        "sinh(1000) is +inf");

    ASSERT_TRUE(&ctx,
        ml_isinf(ml_sinh(-1000.0)) && ml_sinh(-1000.0) < 0.0,
        "sinh(-1000) is -inf");

    ASSERT_TRUE(&ctx,
        ml_isinf(ml_cosh(1000.0)) && ml_cosh(1000.0) > 0.0,
        "cosh(1000) is +inf");

    ASSERT_TRUE(&ctx,
        ml_isinf(ml_cosh(-1000.0)) && ml_cosh(-1000.0) > 0.0,
        "cosh(-1000) is +inf");

    return ml_test_summary(&ctx);
}
