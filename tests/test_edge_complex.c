/* v11S CLOSURE IP-20: edge complex tests */
#include "test_harness.h"
#include "ml_complex.h"
int main (void) {
    ml_test_ctx_t ctx;
    ml_test_init (&ctx, "Edge Complex");

    cplx big = {1e308, 1e308};
    double abs_big = ml_cplx_abs (big);

    ASSERT_TRUE (&ctx, ml_isfinite (abs_big), "complex abs overflow-safe");
    ASSERT_NEAR (&ctx, abs_big, 1.4142135623730951e308, 1e296, "complex abs value");

    cplx zero = {0.0, 0.0};
    ASSERT_NEAR (&ctx, ml_cplx_abs (zero), 0.0, 0.0, "complex abs zero");

    ASSERT_NEAR (&ctx, ml_cplx_arg ((cplx) {1.0, 0.0}), 0.0, 1e-15, "arg(1)");
    ASSERT_NEAR (&ctx, ml_cplx_arg ((cplx) {0.0, 1.0}), ML_PI / 2.0, 1e-15, "arg(i)");
    ASSERT_NEAR (&ctx, ml_cplx_arg ((cplx) {-1.0, 0.0}), ML_PI, 1e-15, "arg(-1)");
    ASSERT_NEAR (&ctx, ml_cplx_arg ((cplx) {0.0, -1.0}), -ML_PI / 2.0, 1e-15, "arg(-i)");

    ASSERT_TRUE (&ctx, ml_isnan (ml_cplx_arg ((cplx) {ml_make_nan (), 0.0})), "arg(NaN) is NaN");

    cplx div_zero = ml_cplx_div ((cplx) {1.0, 0.0}, (cplx) {0.0, 0.0});
    ASSERT_TRUE (&ctx, ml_isnan (div_zero.real) && ml_isnan (div_zero.imag), "complex div by zero gives NaN");

    cplx exp_nan = ml_cplx_exponential ((cplx) {ml_make_nan (), 0.0});
    ASSERT_TRUE (&ctx, ml_isnan (exp_nan.real), "complex exp NaN guard");

    cplx log_zero = ml_cplx_logarithm ((cplx) {0.0, 0.0});
    ASSERT_TRUE (&ctx, ml_isinf (log_zero.real) && log_zero.real < 0.0, "complex log(0) real is -inf");

    /* ---- ml_cplx_sqrt: infinite values (C99 Annex G) ---- */
    {
        /* csqrt(x + i*Inf) = +Inf + i*Inf for any finite x */
        cplx r1 = ml_cplx_sqrt ((cplx) {-5.0, ml_make_inf (0)});
        ASSERT_TRUE (&ctx, ml_isinf (r1.real) && r1.real > 0.0, "csqrt(-5+iInf): real +inf");
        ASSERT_TRUE (&ctx, ml_isinf (r1.imag) && r1.imag > 0.0, "csqrt(-5+iInf): imag +inf");

        cplx r2 = ml_cplx_sqrt ((cplx) {5.0, -ml_make_inf (0)});
        ASSERT_TRUE (&ctx, ml_isinf (r2.real) && r2.real > 0.0, "csqrt(5-iInf): real +inf");
        ASSERT_TRUE (&ctx, ml_isinf (r2.imag) && r2.imag < 0.0, "csqrt(5-iInf): imag -inf");

        /* csqrt(Inf + i*Inf) = +Inf + i*Inf */
        cplx r3 = ml_cplx_sqrt ((cplx) {ml_make_inf (0), ml_make_inf (0)});
        ASSERT_TRUE (&ctx, ml_isinf (r3.real) && r3.real > 0.0, "csqrt(Inf+iInf): real +inf");
        ASSERT_TRUE (&ctx, ml_isinf (r3.imag) && r3.imag > 0.0, "csqrt(Inf+iInf): imag +inf");

        /* csqrt(Inf + i*finite) = +Inf + i*0 with sign of imag */
        cplx r4 = ml_cplx_sqrt ((cplx) {ml_make_inf (0), 3.0});
        ASSERT_TRUE (&ctx, ml_isinf (r4.real) && r4.real > 0.0, "csqrt(Inf+3i): real +inf");
        ASSERT_TRUE (&ctx, r4.imag == 0.0 && ml_signbit (r4.imag) == 0, "csqrt(Inf+3i): imag +0");

        cplx r5 = ml_cplx_sqrt ((cplx) {ml_make_inf (0), -3.0});
        ASSERT_TRUE (&ctx, ml_isinf (r5.real) && r5.real > 0.0, "csqrt(Inf-3i): real +inf");
        ASSERT_TRUE (&ctx, r5.imag == 0.0 && ml_signbit (r5.imag) != 0, "csqrt(Inf-3i): imag -0");

        /* csqrt(-Inf + i*finite) = +0 + i*Inf with sign of imag */
        cplx r6 = ml_cplx_sqrt ((cplx) {-ml_make_inf (0), 3.0});
        ASSERT_TRUE (&ctx, r6.real == 0.0 && ml_signbit (r6.real) == 0, "csqrt(-Inf+3i): real +0");
        ASSERT_TRUE (&ctx, ml_isinf (r6.imag) && r6.imag > 0.0, "csqrt(-Inf+3i): imag +inf");

        cplx r7 = ml_cplx_sqrt ((cplx) {-ml_make_inf (0), -3.0});
        ASSERT_TRUE (&ctx, r7.real == 0.0 && ml_signbit (r7.real) == 0, "csqrt(-Inf-3i): real +0");
        ASSERT_TRUE (&ctx, ml_isinf (r7.imag) && r7.imag < 0.0, "csqrt(-Inf-3i): imag -inf");
    }

    /* ---- ml_cplx_sqrt: zeros and signed zeros ---- */
    {
        /* csqrt(-0,+0) = (+0,+0) */
        cplx r = ml_cplx_sqrt ((cplx) {-0.0, 0.0});
        ASSERT_TRUE (&ctx, r.real == 0.0 && ml_signbit (r.real) == 0, "csqrt(-0,+0): real +0");
        ASSERT_TRUE (&ctx, r.imag == 0.0 && ml_signbit (r.imag) == 0, "csqrt(-0,+0): imag +0");

        /* csqrt(+0,-0) = (+0,-0) */
        r = ml_cplx_sqrt ((cplx) {0.0, -0.0});
        ASSERT_TRUE (&ctx, r.real == 0.0 && ml_signbit (r.real) == 0, "csqrt(+0,-0): real +0");
        ASSERT_TRUE (&ctx, r.imag == 0.0 && ml_signbit (r.imag) != 0, "csqrt(+0,-0): imag -0");

        /* csqrt(-1,+0) = (+0,+1) */
        r = ml_cplx_sqrt ((cplx) {-1.0, 0.0});
        ASSERT_NEAR (&ctx, r.real, 0.0, 1e-15, "csqrt(-1,+0): real ~0");
        ASSERT_NEAR (&ctx, r.imag, 1.0, 1e-15, "csqrt(-1,+0): imag ~1");

        /* csqrt(-1,-0) = (+0,-1) */
        r = ml_cplx_sqrt ((cplx) {-1.0, -0.0});
        ASSERT_NEAR (&ctx, r.real, 0.0, 1e-15, "csqrt(-1,-0): real ~0");
        ASSERT_NEAR (&ctx, r.imag, -1.0, 1e-15, "csqrt(-1,-0): imag ~-1");
    }

    /* ---- ml_cplx_sqrt: branch cuts and quadrants ---- */
    {
        /* sqrt(-4 + 0i) = 0 + 2i */
        cplx r = ml_cplx_sqrt ((cplx) {-4.0, 0.0});
        ASSERT_NEAR (&ctx, r.real, 0.0, 1e-15, "sqrt(-4+0i): real ~0");
        ASSERT_NEAR (&ctx, r.imag, 2.0, 1e-15, "sqrt(-4+0i): imag ~2");

        /* sqrt(-4 - 0i) = 0 - 2i */
        r = ml_cplx_sqrt ((cplx) {-4.0, -0.0});
        ASSERT_NEAR (&ctx, r.real, 0.0, 1e-15, "sqrt(-4-0i): real ~0");
        ASSERT_NEAR (&ctx, r.imag, -2.0, 1e-15, "sqrt(-4-0i): imag ~-2");

        /* sqrt(4 + 0i) = 2 + 0i */
        r = ml_cplx_sqrt ((cplx) {4.0, 0.0});
        ASSERT_NEAR (&ctx, r.real, 2.0, 1e-15, "sqrt(4+0i): real ~2");
        ASSERT_NEAR (&ctx, r.imag, 0.0, 1e-15, "sqrt(4+0i): imag ~0");

        /* Quadrant I: sqrt(3+4i) = 2+i */
        r = ml_cplx_sqrt ((cplx) {3.0, 4.0});
        ASSERT_NEAR (&ctx, r.real, 2.0, 1e-15, "sqrt(3+4i): real ~2");
        ASSERT_NEAR (&ctx, r.imag, 1.0, 1e-15, "sqrt(3+4i): imag ~1");

        /* Verify round-trip: sqrt(z)^2 == z */
        cplx z = {0.5, -1.25};
        r = ml_cplx_sqrt (z);
        cplx sq = ml_cplx_mul (r, r);
        ASSERT_NEAR (&ctx, sq.real, z.real, 1e-14, "sqrt round-trip: real");
        ASSERT_NEAR (&ctx, sq.imag, z.imag, 1e-14, "sqrt round-trip: imag");
    }

    /* ---- ml_cplx_power ---- */
    {
        /* z^0 = 1 for any z (including 0) */
        cplx r = ml_cplx_power ((cplx) {2.0, 3.0}, (cplx) {0.0, 0.0});
        ASSERT_NEAR (&ctx, r.real, 1.0, 1e-15, "z^0: real ~1");
        ASSERT_NEAR (&ctx, r.imag, 0.0, 1e-15, "z^0: imag ~0");

        /* z^1 = z */
        r = ml_cplx_power ((cplx) {2.0, 3.0}, (cplx) {1.0, 0.0});
        ASSERT_NEAR (&ctx, r.real, 2.0, 1e-15, "z^1: real ~2");
        ASSERT_NEAR (&ctx, r.imag, 3.0, 1e-15, "z^1: imag ~3");

        /* 0^positive = 0 */
        r = ml_cplx_power ((cplx) {0.0, 0.0}, (cplx) {2.0, 0.0});
        ASSERT_NEAR (&ctx, r.real, 0.0, 1e-15, "0^2: real ~0");
        ASSERT_NEAR (&ctx, r.imag, 0.0, 1e-15, "0^2: imag ~0");

        /* i^2 = -1 */
        r = ml_cplx_power ((cplx) {0.0, 1.0}, (cplx) {2.0, 0.0});
        ASSERT_NEAR (&ctx, r.real, -1.0, 1e-14, "i^2: real ~-1");
        ASSERT_NEAR (&ctx, r.imag, 0.0, 1e-14, "i^2: imag ~0");

        /* (1+i)^2 = 2i */
        r = ml_cplx_power ((cplx) {1.0, 1.0}, (cplx) {2.0, 0.0});
        ASSERT_NEAR (&ctx, r.real, 0.0, 1e-14, "(1+i)^2: real ~0");
        ASSERT_NEAR (&ctx, r.imag, 2.0, 1e-14, "(1+i)^2: imag ~2");
    }

    /* ---- ml_cplx_nth_root ---- */
    {
        /* 4th roots of 1: 1, i, -1, -i */
        cplx roots[4];
        for (int k = 0; k < 4; k++) {
            roots[k] = ml_cplx_nth_root ((cplx) {1.0, 0.0}, 4, k);
        }
        /* Root 0: 1 */
        ASSERT_NEAR (&ctx, roots[0].real, 1.0, 1e-14, "4th root of 1, k=0: real ~1");
        ASSERT_NEAR (&ctx, roots[0].imag, 0.0, 1e-14, "4th root of 1, k=0: imag ~0");
        /* Root 1: i */
        ASSERT_NEAR (&ctx, roots[1].real, 0.0, 1e-14, "4th root of 1, k=1: real ~0");
        ASSERT_NEAR (&ctx, roots[1].imag, 1.0, 1e-14, "4th root of 1, k=1: imag ~1");
        /* Root 2: -1 */
        ASSERT_NEAR (&ctx, roots[2].real, -1.0, 1e-14, "4th root of 1, k=2: real ~-1");
        ASSERT_NEAR (&ctx, roots[2].imag, 0.0, 1e-14, "4th root of 1, k=2: imag ~0");
        /* Root 3: -i */
        ASSERT_NEAR (&ctx, roots[3].real, 0.0, 1e-14, "4th root of 1, k=3: real ~0");
        ASSERT_NEAR (&ctx, roots[3].imag, -1.0, 1e-14, "4th root of 1, k=3: imag ~-1");

        /* Verify round-trip: (nth_root(z,n,k))^n == z */
        cplx z = {-8.0, 0.0};
        cplx root = ml_cplx_nth_root (z, 3, 0);
        cplx cube = ml_cplx_mul (ml_cplx_mul (root, root), root);
        ASSERT_NEAR (&ctx, cube.real, z.real, 1e-12, "cube root round-trip: real");
        ASSERT_NEAR (&ctx, cube.imag, z.imag, 1e-12, "cube root round-trip: imag");
    }

    /* ---- ml_cplx_sqrt: NaN propagation ---- */
    {
        cplx r = ml_cplx_sqrt ((cplx) {ml_make_nan (), 0.0});
        ASSERT_TRUE (&ctx, ml_isnan (r.real) || ml_isnan (r.imag), "csqrt(NaN,0): NaN propagated");

        r = ml_cplx_sqrt ((cplx) {0.0, ml_make_nan ()});
        ASSERT_TRUE (&ctx, ml_isnan (r.real) || ml_isnan (r.imag), "csqrt(0,NaN): NaN propagated");
    }

    return ml_test_summary (&ctx);
}
