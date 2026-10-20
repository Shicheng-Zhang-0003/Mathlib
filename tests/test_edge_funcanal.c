/* Doctorate funcanal pins: Banach, Neumann, Hilbert, Rayleigh. */
#include "test_harness.h"
#include "ml_funcanal.h"
#include "ml_core.h"
static double t_cos (double x) { return 0.5 * x + 0.25; }
int main (void) {
    ml_test_ctx_t ctx;
    ml_test_init (&ctx, "Edge Funcanal");

    /* Banach: x = x/2 + 1/4 -> x* = 1/2, L = 1/2 certified. */
    {
        double xs = 0, bd = 0;
        ml_status_t st = ml_banach_fixed_point (t_cos, 0.0, 0.5, 1e-12, 10000, &xs, &bd);
        ASSERT_TRUE (&ctx, st == ML_SUCCESS, "banach status");
        ASSERT_NEAR (&ctx, xs, 0.5, 1e-9, "banach fixed point");
        ASSERT_TRUE (&ctx, bd >= 0.0 && bd < 1e-6, "banach bound sane");
        ASSERT_TRUE (&ctx, ml_banach_fixed_point (t_cos, 0.0, 1.0, 1e-12, 100, &xs, NULL) == ML_ERR_INVALID_ARG, "banach L=1");
        ASSERT_TRUE (&ctx, ml_banach_fixed_point (NULL, 0.0, 0.5, 1e-12, 100, &xs, NULL) == ML_ERR_INVALID_ARG, "banach null");
    }
    /* Neumann: A = [[0.8, 0.1],[0.1, 0.8]] -> inv = [[1.316, -0.164],...]. */
    {
        double A [] = {0.8, 0.1, 0.1, 0.8};
        double out [4] = {0, 0, 0, 0};
        ml_status_t st = ml_neumann_inverse (A, 2, 200, out);
        ASSERT_TRUE (&ctx, st == ML_SUCCESS, "neumann status");
        ASSERT_NEAR (&ctx, out [0], 1.2698412698, 1e-6, "neumann 00");
        ASSERT_NEAR (&ctx, out [1], -0.1587301587, 1e-6, "neumann 01");
        ASSERT_TRUE (&ctx, ml_neumann_inverse ((double []) {2, 0, 0, 2}, 2, 50, out) == ML_ERR_INVALID_ARG, "neumann divergent");
    }
    /* Hilbert: project (1,1,1) onto span{(1,0,0),(0,1,0)} -> (1,1,0). */
    {
        double V [] = {1, 0, 0, 1, 0, 0};
        double y [] = {1, 1, 1};
        double c [2] = {0, 0};
        ml_status_t st = ml_hilbert_project (V, y, 3, 2, c);
        ASSERT_TRUE (&ctx, st == ML_SUCCESS, "hilbert status");
        ASSERT_NEAR (&ctx, c [0], 1.0, 1e-9, "hilbert c0");
        ASSERT_NEAR (&ctx, c [1], 1.0, 1e-9, "hilbert c1");
        ASSERT_TRUE (&ctx, ml_hilbert_project (V, V, 3, 2, c) == ML_SUCCESS, "hilbert self");
    }
    /* Rank loss: collinear columns -> SINGULAR. */
    {
        double V [] = {1, 2, 2, 4, 3, 6};
        double y [] = {1, 0, 0};
        double c [2] = {0, 0};
        ASSERT_TRUE (&ctx, ml_hilbert_project (V, y, 3, 2, c) == ML_ERR_SINGULAR, "hilbert rank loss");
    }
    /* Rayleigh: diag(3,1) -> 3. */
    {
        double A [] = {3, 0, 0, 1};
        double lam = 0;
        ml_status_t st = ml_rayleigh_max (A, 2, &lam);
        ASSERT_TRUE (&ctx, st == ML_SUCCESS, "rayleigh status");
        ASSERT_NEAR (&ctx, lam, 3.0, 1e-9, "rayleigh max");
    }

    return ml_test_summary (&ctx);
}
