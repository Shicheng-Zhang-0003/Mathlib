/* Despot Round-4 regression pins: zeta(0), Hurwitz EM, complex LD paths,
 * remainder big-quotient, gamma P/Q at Inf.
 * Ground truth: mpmath 80-dps (see DESPOT_R4_AUDIT.md); tolerances are
 * regression guards, not oracle certification (PRECISION_CONTRACT.md). */
#include "test_harness.h"
#include "ml_complex.h"
#include "ml_integral.h"
#include "ml_analytic_nt.h"
#include "ml_core.h"
int main (void) {
    ml_test_ctx_t ctx;
    ml_test_init (&ctx, "Edge Despot R4");

    /* F1: zeta(0) = -1/2 (was 0 via trivial-zero guard including s==0). */
    ASSERT_NEAR (&ctx, ml_zeta (0.0), -0.5, 0.0, "zeta(0)=-1/2");
    ASSERT_NEAR (&ctx, ml_zeta (-0.0), -0.5, 0.0, "zeta(-0)=-1/2");
    ASSERT_NEAR (&ctx, ml_zeta (-2.0), 0.0, 0.0, "zeta(-2) trivial zero");
    ASSERT_NEAR (&ctx, ml_zeta (-4.0), 0.0, 0.0, "zeta(-4) trivial zero");
    ASSERT_NEAR (&ctx, ml_zeta (2.0), 1.6449340668482264, 1e-15, "zeta(2)=pi^2/6");

    /* Hurwitz EM N=200 B2..B12: 0 ULP vs mpmath 80-dps. */
    ASSERT_NEAR (&ctx, ml_hurwitz_zeta (1.1, 1.0).real, 10.584448464950801, 1e-12, "hurwitz(1.1,1)");
    ASSERT_NEAR (&ctx, ml_hurwitz_zeta (1.5, 2.0).real, 1.6123753486854884, 1e-14, "hurwitz(1.5,2)");
    ASSERT_NEAR (&ctx, ml_hurwitz_zeta (0.9, 2.0).real, -10.430114019402255, 1e-11, "hurwitz(0.9,2) continuation");
    ASSERT_TRUE (&ctx, ml_isnan (ml_hurwitz_zeta (-1.0, 2.0).real), "hurwitz(-1,2) fail-loud NaN stub");

    /* Complex log: 0.5*logl(x^2+y^2) LD path (was 88x amplification). */
    {
        cplx l = ml_cplx_logarithm ((cplx) {1.0, 0.1});
        ASSERT_NEAR (&ctx, l.real, 0.0049751654265840417, 1e-17, "clog(1+0.1i).real");
        ASSERT_NEAR (&ctx, l.imag, 0.099668652491162038, 1e-15, "clog(1+0.1i).imag");
    }

    /* Complex sqrt Hull-Beebe division form (was im=0 for 1e300+1e200i). */
    {
        cplx s = ml_cplx_sqrt ((cplx) {1e300, 1e200});
        ASSERT_NEAR (&ctx, s.real, 1e150, 1e135, "csqrt(1e300+1e200i).real");
        ASSERT_TRUE (&ctx, s.imag > 4e49 && s.imag < 6e49, "csqrt(1e300+1e200i).imag ~5e49");
    }

    /* Complex pow LD w (was 885 ULP): (1+0.1i)^37 vs mpmath. */
    {
        cplx p = ml_cplx_power ((cplx) {1.0, 0.1}, (cplx) {37.0, 0.0});
        ASSERT_NEAR (&ctx, p.real, -1.0272441314090439, 1e-14, "cpow real");
        ASSERT_NEAR (&ctx, p.imag, -0.62437646169052141, 1e-14, "cpow imag");
    }
    /* pow(0,0)=1 edge still holds with LD guard present. */
    {
        cplx p = ml_cplx_power ((cplx) {0.0, 0.0}, (cplx) {0.0, 0.0});
        ASSERT_NEAR (&ctx, p.real, 1.0, 0.0, "cpow(0,0)=1 real");
        ASSERT_NEAR (&ctx, p.imag, 0.0, 0.0, "cpow(0,0)=1 imag");
    }

    /* Remainder big-quotient (was 3.5e283 vs 160). */
    ASSERT_NEAR (&ctx, ml_remainder (1e300, 1000.0), remainder (1e300, 1000.0), 0.0, "remainder(1e300,1000) bit-exact");
    ASSERT_NEAR (&ctx, ml_remainder (5.3, 2.0), remainder (5.3, 2.0), 0.0, "remainder(5.3,2)");

    /* gamma P/Q at +Inf (was NaN via erf(1e300)). */
    ASSERT_NEAR (&ctx, ml_gamma_p (0.5, 1e300), 1.0, 0.0, "gamma_p(1/2,inf)=1");
    ASSERT_NEAR (&ctx, ml_gamma_q (0.5, 1e300), 0.0, 0.0, "gamma_q(1/2,inf)=0");

    return ml_test_summary (&ctx);
}
