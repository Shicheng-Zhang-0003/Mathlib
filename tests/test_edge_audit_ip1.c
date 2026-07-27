/* v11S AUDIT IP-1: regression tests */
#include "test_harness.h"
#include "ml_core.h"
#include "ml_trig.h"
#include "ml_complex.h"
#include "ieee754.h"
#include "fast_math.h"
#include "simd_batch.h"
#include <stdint.h>
#include <string.h>

int main(void) {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Edge Audit IP-1");

    /* atan2 NaN / signed-zero / infinity behavior */
    ASSERT_TRUE(&ctx, ml_isnan(ml_atan2(0.0, ml_make_nan())), "atan2(0,NaN) is NaN");
    ASSERT_TRUE(&ctx, ml_isnan(ml_atan2(ml_make_nan(), 0.0)), "atan2(NaN,0) is NaN");
    ASSERT_TRUE(&ctx, ml_signbit(ml_atan2(-0.0, ml_make_inf(0))) != 0, "atan2(-0,+inf) preserves sign");
    ASSERT_TRUE(&ctx, ml_signbit(ml_atan2(-1.0, ml_make_inf(0))) != 0, "atan2(-finite,+inf) preserves sign");
    ASSERT_NEAR(&ctx, ml_atan2(0.0, -1.0), ML_PI, 1e-15, "atan2(+0,-1)");
    ASSERT_NEAR(&ctx, ml_atan2(-0.0, -1.0), -ML_PI, 1e-15, "atan2(-0,-1)");
    ASSERT_TRUE(&ctx, ml_signbit(ml_atan2(-0.0, 1.0)) != 0, "atan2(-0,+1) preserves sign");

    /* complex power special cases */
    cplx zero = {0.0, 0.0};
    cplx one = {1.0, 0.0};
    cplx nanv = {ml_make_nan(), 0.0};

    cplx p00 = ml_cplx_power(zero, zero);
    ASSERT_TRUE(&ctx, p00.real == 1.0 && p00.imag == 0.0, "cplx pow(0,0) == 1");

    cplx p1n = ml_cplx_power(one, nanv);
    ASSERT_TRUE(&ctx, p1n.real == 1.0 && p1n.imag == 0.0, "cplx pow(1,NaN) == 1");

    /* experimental ieee754 helpers are now safe */
    ASSERT_TRUE(&ctx, ml_isnan(logarithm_ieee754(-1.0)), "ieee log(-1) is NaN");
    ASSERT_TRUE(&ctx, ml_isinf(logarithm_ieee754(0.0)) && logarithm_ieee754(0.0) < 0.0, "ieee log(0) is -inf");
    ASSERT_NEAR(&ctx, logarithm_ieee754(1.0), 0.0, 1e-15, "ieee log(1)");
    ASSERT_NEAR(&ctx, logarithm_ieee754(2.0), ML_LN2, 1e-12, "ieee log(2)");
    ASSERT_NEAR(&ctx, logarithm_ieee754(0.5), -ML_LN2, 1e-12, "ieee log(0.5)");

    ASSERT_NEAR(&ctx, exponential_ieee754(1.0), ML_E, 1e-10, "ieee exp(1)");
    ASSERT_TRUE(&ctx, ml_isinf(exponential_ieee754(1000.0)) && exponential_ieee754(1000.0) > 0.0, "ieee exp(1000) is +inf");
    ASSERT_TRUE(&ctx, exponential_ieee754(-1000.0) == 0.0, "ieee exp(-1000) is 0");

    /* fast math subnormal fallback */
    uint64_t min_bits = 1;
    double min_sub;
    memcpy(&min_sub, &min_bits, sizeof(double));
    ASSERT_TRUE(&ctx, ml_isfinite(ml_fast_rsqrt(min_sub)), "fast rsqrt subnormal finite");

    /* SIMD batch NaN propagation is now unconditional after P1 guard */
    double in[4];
    double out[4];

    in[0] = 1.0;
    in[1] = 4.0;
    in[2] = 0.0;
    in[3] = ml_make_nan();

    ml_simd_batch_rsqrt(in, out);

    /* MATHLIB_CLOSURE_P1_SIMD_BATCH_GUARD_TEST */
    ASSERT_TRUE(&ctx, ml_isnan(out[3]), "batch rsqrt propagates NaN");

    /* fmod exactness basics */
    ASSERT_NEAR(&ctx, ml_fmod(5.0, 3.0), 2.0, 1e-15, "fmod(5,3)");
    ASSERT_NEAR(&ctx, ml_fmod(-5.0, 3.0), -2.0, 1e-15, "fmod(-5,3)");
    ASSERT_TRUE(&ctx, ml_isnan(ml_fmod(1.0, 0.0)), "fmod(x,0) is NaN");
    ASSERT_TRUE(&ctx, ml_signbit(ml_fmod(-0.0, 1.0)) != 0, "fmod(-0,1) preserves sign");

    return ml_test_summary(&ctx);
}
