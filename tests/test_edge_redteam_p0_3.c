/* MATHLIB_CLOSURE_P2_P0_3_TEST */
/* Red Team P2 P0-3: ml_exp threshold regression tests */
#include "test_harness.h"
#include "ml_core.h"
#include "ml_exp_log.h"

int main(void) {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Edge Red Team P0-3");

    /*
     * Old overflow threshold was too low:
     *
     *   if (x > 709.78) return inf;
     *
     * The true approximate overflow boundary is:
     *
     *   log(DBL_MAX) ~= 709.782712893384
     *
     * Therefore values slightly above 709.78 but below log(DBL_MAX)
     * must remain finite.
     */
    ASSERT_TRUE(&ctx,
        ml_isfinite(ml_exp(709.781)),
        "exp(709.781) must remain finite");

    ASSERT_TRUE(&ctx,
        ml_isfinite(ml_exp(709.782)),
        "exp(709.782) must remain finite");

    /*
     * Values above log(DBL_MAX) must overflow.
     */
    ASSERT_TRUE(&ctx,
        ml_isinf(ml_exp(709.783)) && ml_exp(709.783) > 0.0,
        "exp(709.783) must overflow to +inf");

    ASSERT_TRUE(&ctx,
        ml_isinf(ml_exp(710.0)) && ml_exp(710.0) > 0.0,
        "exp(710) must overflow to +inf");

    /*
     * Old underflow threshold was too high:
     *
     *   if (x < -745.13) return 0;
     *
     * The approximate underflow-to-zero boundary is closer to:
     *
     *   -745.133219101941
     *
     * Therefore -745.132 should still produce a positive subnormal
     * or normal result, not zero.
     */
    ASSERT_TRUE(&ctx,
        ml_exp(-745.132) > 0.0,
        "exp(-745.132) must not underflow prematurely");

    ASSERT_TRUE(&ctx,
        ml_exp(-745.2) == 0.0,
        "exp(-745.2) underflows to zero");

    /*
     * Basic special-case sanity must remain intact.
     */
    ASSERT_TRUE(&ctx,
        ml_exp(0.0) == 1.0,
        "exp(0) == 1");

    ASSERT_TRUE(&ctx,
        ml_isnan(ml_exp(ml_make_nan())),
        "exp(NaN) is NaN");

    ASSERT_TRUE(&ctx,
        ml_isinf(ml_exp(ml_make_inf(0))) && ml_exp(ml_make_inf(0)) > 0.0,
        "exp(+inf) is +inf");

    ASSERT_TRUE(&ctx,
        ml_exp(-ml_make_inf(0)) == 0.0,
        "exp(-inf) is 0");

    return ml_test_summary(&ctx);
}
