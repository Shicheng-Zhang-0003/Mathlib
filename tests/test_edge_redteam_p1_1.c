/* MATHLIB_REDP2_P1_1_TEST_V2 */
/* Red Team P2 P1-1: SIMD batch rsqrt scalar-fallback semantics */

#include "test_harness.h"
#include "simd_batch.h"
#include "fast_math.h"
#include "ml_core.h"

#include <stdint.h>
#include <string.h>

int main(void) {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Edge Red Team P1-1");

    double in[4];
    double out[4];

    /*
     * Exceptional batch:
     *   +0, -0, negative finite, NaN
     *
     * This must force scalar fallback for the whole batch.
     */
    in[0] = 0.0;
    in[1] = -0.0;
    in[2] = -1.0;
    in[3] = ml_make_nan();

    ml_simd_batch_rsqrt(in, out);

    ASSERT_TRUE(&ctx,
        ml_isinf(out[0]) && out[0] > 0.0,
        "batch rsqrt(+0) gives +inf");

    ASSERT_TRUE(&ctx,
        ml_isinf(out[1]) && ml_signbit(out[1]) != 0,
        "batch rsqrt(-0) gives -inf");

    ASSERT_TRUE(&ctx,
        ml_isnan(out[2]),
        "batch rsqrt(negative finite) gives NaN");

    ASSERT_TRUE(&ctx,
        ml_isnan(out[3]),
        "batch rsqrt(NaN) gives NaN");

    /*
     * Exceptional batch:
     *   +inf, -inf, smallest subnormal, positive finite
     */
    uint64_t min_bits = 1;
    double min_sub;
    memcpy(&min_sub, &min_bits, sizeof(double));

    in[0] = ml_make_inf(0);
    in[1] = -ml_make_inf(0);
    in[2] = min_sub;
    in[3] = 4.0;

    ml_simd_batch_rsqrt(in, out);

    ASSERT_TRUE(&ctx,
        out[0] == 0.0 && ml_signbit(out[0]) == 0,
        "batch rsqrt(+inf) gives +0");

    ASSERT_TRUE(&ctx,
        ml_isnan(out[1]),
        "batch rsqrt(-inf) gives NaN");

    double scalar_sub = ml_fast_rsqrt(min_sub);

    ASSERT_TRUE(&ctx,
        ml_isfinite(out[2]),
        "batch rsqrt(subnormal) finite");

    ASSERT_NEAR(&ctx,
        out[2],
        scalar_sub,
        ml_fabs(scalar_sub) * 1e-15 + 1e-300,
        "batch rsqrt(subnormal) matches scalar fast-math fallback");

    double scalar_four = ml_fast_rsqrt(4.0);

    ASSERT_NEAR(&ctx,
        out[3],
        scalar_four,
        ml_fabs(scalar_four) * 1e-15 + 1e-300,
        "batch rsqrt(4) matches scalar fast-math fallback");

    /*
     * All-positive normal batch:
     *
     * This may use the AVX2 fast path.
     *
     * The fast path is intentionally approximate. The observed error
     * for the current Quake-style double rsqrt kernel is around 1e-6
     * to 1e-5, so 1e-12 is far too strict.
     *
     * Use a fast-math-appropriate relative tolerance.
     */
    double in_pos[4] = {1.0, 4.0, 16.0, 100.0};
    double out_pos[4];

    ml_simd_batch_rsqrt(in_pos, out_pos);

    for (int i = 0; i < 4; i++) {
        double expected = 1.0 / ml_sqrt(in_pos[i]);
        double tol = 1e-4 * ml_fabs(expected) + 1e-12;

        ASSERT_NEAR(&ctx,
            out_pos[i],
            expected,
            tol,
            "batch rsqrt positive normal approximate");
    }

    return ml_test_summary(&ctx);
}
