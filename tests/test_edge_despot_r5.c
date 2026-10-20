/* Despot Round-5 regression pins: every R5 programming/operational fix.
 * Ground truth: mpmath 80-dps / DLMF identities / closed-form physics
 * (see docs/DESPOT_R5_AUDIT.md); tolerances are regression guards, not
 * oracle certification (see PRECISION_CONTRACT.md rule 3). */
#include "test_harness.h"
#include "ml_core.h"
#include "ml_combinatorics.h"
#include "ml_numbertheory.h"
#include "ml_numerical.h"
#include "ml_linalg.h"
#include "ml_tensor.h"
#include "ml_spectral.h"
#include "ml_ode.h"
#include "ml_ode_sys.h"
#include "ml_sde.h"
#include "ml_mcmc.h"
#include "ml_info.h"
#include "ml_transforms.h"
#include "ml_analytic_nt.h"
#include "ml_integral.h"
#include "ml_fixed_point.h"
#include "cpu_dispatch.h"
#include "fft.h"
static double f_exp (double t, double y) {
    (void) t;
    return y;
}
static void sys_exp (double t, const double *y, double *dydt, int n, void *ctx) {
    (void) t;
    (void) ctx;
    for (int i = 0; i < n; i++) dydt [i] = y [i];
}
static double drift_c (double t, double x, void *c) {
    (void) t;
    (void) x;
    (void) c;
    return 0.0;
}
static double diff_c (double t, double x, void *c) {
    (void) t;
    (void) x;
    (void) c;
    return 1.0;
}
int main (void) {
    ml_test_ctx_t ctx;
    ml_test_init (&ctx, "Edge Despot R5");

    /* A1: cpu dispatch never SIGILLs; scalar fallback is safe. */
    {
        int avx = ml_cpu_has_avx2 ();
        ASSERT_TRUE (&ctx, avx == 0 || avx == 1, "cpu_has_avx2 boolean");
        double A [] = {1, 0, 0, 1}, B [] = {5, 6, 7, 8}, C [] = {0, 0, 0, 0};
        ml_matmul (A, B, C, 2);
        ASSERT_NEAR (&ctx, C [0], 5.0, 0.0, "matmul scalar fallback");
        ASSERT_NEAR (&ctx, C [3], 8.0, 0.0, "matmul identity");
    }

    /* A2: cross3 alias-safe (was wrong y,z for out==a). */
    {
        double a [] = {1, 0, 0}, b [] = {0, 1, 0};
        ml_cross3 (a, b, a);
        ASSERT_NEAR (&ctx, a [0], 0.0, 0.0, "cross3 alias x");
        ASSERT_NEAR (&ctx, a [1], 0.0, 0.0, "cross3 alias y");
        ASSERT_NEAR (&ctx, a [2], 1.0, 0.0, "cross3 alias z");
    }

    /* A3: transpose in-place square (was garbage without temp). */
    {
        double d [] = {1, 2, 3, 4};
        ml_tensor_view_t A = {d, 2, 2}, T = {d, 2, 2};
        ml_transpose (A, T);
        ASSERT_NEAR (&ctx, d [1], 3.0, 0.0, "transpose in-place (0,1)");
        ASSERT_NEAR (&ctx, d [2], 2.0, 0.0, "transpose in-place (1,0)");
    }

    /* A4: egcd iterative, overflow-checked (was recursive UB). */
    {
        int64_t x = 0, y = 0;
        int64_t g = ml_egcd (48, 18, &x, &y);
        ASSERT_TRUE (&ctx, g == 6, "egcd(48,18)=6");
        ASSERT_TRUE (&ctx, 48 * x + 18 * y == 6, "egcd bezout");
        g = ml_egcd (-5, 0, &x, &y);
        ASSERT_TRUE (&ctx, g == 5 && x == -1 && y == 0, "egcd(-5,0) sign");
        g = ml_egcd (INT64_MIN, 1, &x, &y);
        ASSERT_TRUE (&ctx, g == -1, "egcd(INT64_MIN) sentinel");
    }

    /* A5: fft guards (1<<24 signed edge -> 16777216, boolean contract). */
    {
        ASSERT_TRUE (&ctx, ml_fft_is_supported (1024) == 1, "fft pow2 supported");
        ASSERT_TRUE (&ctx, ml_fft_is_supported (1000) == 0, "fft non-pow2 refused");
        ASSERT_TRUE (&ctx, ml_fft_is_supported (-4) == 0, "fft negative refused");
    }

    /* B1: status-code precision (CONVERGENCE/OVERFLOW exist and are used). */
    {
        ASSERT_TRUE (&ctx, ML_ERR_CONVERGENCE == -6, "CONVERGENCE code");
        ASSERT_TRUE (&ctx, ML_ERR_OVERFLOW == -5, "OVERFLOW code");
        /* True singularity keeps SINGULAR (ABI preserved). */
        double Az [] = {0, 0, 0, 0}, bz [] = {0, 0}, xz [] = {0, 0};
        uint8_t wbuf [4096];
        ml_workspace_t ws;
        ml_workspace_init (&ws, wbuf, sizeof (wbuf));
        ml_tensor_view_t Av = {Az, 2, 2};
        ASSERT_TRUE (&ctx, ml_solve (Av, bz, xz, &ws) == ML_ERR_SINGULAR, "zero matrix still SINGULAR");
    }

    /* B2: DoS caps fail loud, not hang/OOM. */
    {
        ASSERT_TRUE (&ctx, ml_majorizes ((double []) {1}, (double []) {1}, 1048577, 1e-12) == 0, "majorizes >1M cap");
        {
            static double big [8];
            for (int i = 0; i < 8; i++) big [i] = (double) i;
            ASSERT_TRUE (&ctx, ml_isnan (ml_ess (big, 100001)), "ess >100k NaN cap");
        }
        ASSERT_TRUE (&ctx, ml_isnan (ml_parseval_energy ((double []) {1}, 1000001)), "parseval >1M NaN cap");
        /* DCT/DST void caps: must not crash on oversize. */
        {
            static double xin [8], xout [8];
            for (int i = 0; i < 8; i++) xin [i] = 1.0;
            ml_dct2 (xin, xout, 5000);
            ml_dct3 (xin, xout, 5000);
            ml_dst2 (xin, xout, 5000);
            ASSERT_TRUE (&ctx, 1, "dct/dst oversize no-crash");
        }
    }

    /* B3: RNG warmup deterministic (same seed => same stream). */
    {
        double s1 = ml_brownian_bridge_sample (0.0, 1.0, 1.0, 0.3, 42ULL);
        double s2 = ml_brownian_bridge_sample (0.0, 1.0, 1.0, 0.3, 42ULL);
        ASSERT_NEAR (&ctx, s1, s2, 0.0, "bridge_sample deterministic");
        double e1 = 0, e2 = 0;
        ml_sde_euler_maruyama_ctx (drift_c, diff_c, 0.0, 0.0, 1.0, 10, 7ULL, &e1, 0);
        ml_sde_euler_maruyama_ctx (drift_c, diff_c, 0.0, 0.0, 1.0, 10, 7ULL, &e2, 0);
        ASSERT_NEAR (&ctx, e1, e2, 0.0, "sde ctx deterministic");
    }

    /* B4: ODE arrival (was t==t1 exact-float stall) + physics truth. */
    {
        double y1 [1] = {0};
        ml_status_t st = ml_ode_dp5_sys (sys_exp, 0.0, (double []) {1.0}, 1, 1.0, 1e-10, 1e-12, y1, 0);
        ASSERT_TRUE (&ctx, st == ML_SUCCESS, "dp5_sys converges");
        ASSERT_NEAR (&ctx, y1 [0], 2.7182818284590451, 1e-8, "dp5_sys exp(1) physics");
        ASSERT_NEAR (&ctx, ml_ode_rk4 (f_exp, 0.0, 1.0, 0.01, 100), 2.7182818284590451, 1e-6, "rk4 exp(1)");
        ASSERT_NEAR (&ctx, ml_suvat_s (0, 9.81, 2.0), 19.62, 1e-12, "SUVAT s=1/2at^2");
        ASSERT_NEAR (&ctx, ml_projectile_range (100.0, 0.7853981633974483, 9.81), 1019.3679918450561, 1e-9, "projectile v0^2/g");
        ASSERT_NEAR (&ctx, ml_shm_x (1.0, 1.0, 6.283185307179586, 0.0), 1.0, 1e-12, "SHM cos period");
    }

    /* B5: Hurwitz large-guard + continuation (was dead clause / NaN stub). */
    {
        ASSERT_NEAR (&ctx, ml_hurwitz_zeta (2.0, 1.0).real, 1.6449340668482264, 1e-12, "hurwitz(2,1)=pi^2/6");
        ASSERT_TRUE (&ctx, ml_isnan (ml_hurwitz_zeta (-1.0, 2.0).real), "hurwitz(-1,2) stub NaN");
        ASSERT_TRUE (&ctx, ml_isfinite (ml_zeta (2.0)), "zeta(2) finite");
    }

    /* B6: combinatorics exactness (pollard cap keeps worst-case bounded). */
    {
        ASSERT_TRUE (&ctx, ml_is_prime64 (18446744073709551557ULL) == 1, "prime64 large prime");
        ASSERT_TRUE (&ctx, ml_gcd (48, 18) == 6, "gcd");
        ASSERT_TRUE (&ctx, ml_phi (12) == 4, "phi(12)=4");
    }

    /* B7: info/entropy physics (fair coin = ln2). */
    {
        double p [] = {0.5, 0.5};
        ASSERT_NEAR (&ctx, ml_entropy (p, 2), 0.69314718055994529, 1e-15, "entropy fair coin");
        ASSERT_NEAR (&ctx, ml_logistic (0.0), 0.5, 0.0, "logistic(0)");
        ASSERT_NEAR (&ctx, ml_softplus (0.0), 0.69314718055994529, 1e-15, "softplus(0)=ln2");
    }

    /* B8: spectral still solves (CONVERGENCE only on true exhaustion). */
    {
        double Ad [] = {4, 0, 0, 9}, b [] = {8, 18}, x [] = {0, 0};
        ml_tensor_view_t A = {Ad, 2, 2};
        ml_status_t st = ml_cg_solve (A, b, x, 100, 1e-12);
        ASSERT_TRUE (&ctx, st == ML_SUCCESS, "cg solves diag");
        ASSERT_NEAR (&ctx, x [0], 2.0, 1e-9, "cg x0");
        ASSERT_NEAR (&ctx, x [1], 2.0, 1e-9, "cg x1");
    }

    return ml_test_summary (&ctx);
}
