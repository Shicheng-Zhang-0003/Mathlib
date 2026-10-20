/* Doctorate diffgeom pins: circle/helix/sphere against closed forms. */
#include "test_harness.h"
#include "ml_diffgeom.h"
#include "ml_core.h"
#include <math.h>
static double circx (double t) { return 2.0 * cos (t); }
static double circy (double t) { return 2.0 * sin (t); }
static double helx (double t) { return cos (t); }
static double hely (double t) { return sin (t); }
static double helz (double t) { return t; }
static double sq (double x) { return x * x; }
static double spx (double u, double v) { return 3.0 * sin (u) * cos (v); }
static double spy (double u, double v) { return 3.0 * sin (u) * sin (v); }
static double spz (double u, double v) {
    (void) v;
    return 3.0 * cos (u);
}
static double plx (double u, double v) {
    (void) v;
    return u;
}
static double ply (double u, double v) {
    (void) u;
    return v;
}
static double plz (double u, double v) {
    (void) u;
    (void) v;
    return 0.0;
}
int main (void) {
    ml_test_ctx_t ctx;
    ml_test_init (&ctx, "Edge Diffgeom");

    ASSERT_NEAR (&ctx, ml_curvature_param (circx, circy, 0.3), 0.5, 1e-9, "circle kappa=1/2");
    ASSERT_NEAR (&ctx, ml_curvature_param (circx, circy, 2.1), 0.5, 1e-9, "circle elsewhere");
    ASSERT_NEAR (&ctx, ml_curvature_graph (sq, 0.0), 2.0, 1e-9, "parabola vertex");
    ASSERT_TRUE (&ctx, ml_isnan (ml_curvature_param (NULL, circy, 0.0)), "null callback NaN");

    /* Helix radius 1 pitch 1: kappa = tau = 1/2. */
    {
        double T [3] = {0, 0, 0}, N [3] = {0, 0, 0}, B [3] = {0, 0, 0};
        double k = 0, t = 0;
        ml_status_t st = ml_frenet_frame (helx, hely, helz, 1.0, T, N, B, &k, &t);
        ASSERT_TRUE (&ctx, st == ML_SUCCESS, "helix status");
        ASSERT_NEAR (&ctx, k, 0.5, 1e-6, "helix kappa");
        ASSERT_NEAR (&ctx, t, 0.5, 1e-3, "helix tau");
        {
            double dot = T [0] * N [0] + T [1] * N [1] + T [2] * N [2];
            double nt = sqrt (N [0] * N [0] + N [1] * N [1] + N [2] * N [2]);
            ASSERT_NEAR (&ctx, dot, 0.0, 1e-9, "T perp N");
            ASSERT_NEAR (&ctx, nt, 1.0, 1e-9, "N unit");
        }
    }
    /* Straight line: SINGULAR, kappa 0, T defined. */
    {
        double T [3] = {0, 0, 0}, N [3] = {0, 0, 0}, B [3] = {0, 0, 0};
        double k = 9, t = 9;
        ml_status_t st = ml_frenet_frame (helx, helx, helx, 0.5, T, N, B, &k, &t);
        ASSERT_TRUE (&ctx, st == ML_ERR_SINGULAR, "straight singular");
        ASSERT_NEAR (&ctx, k, 0.0, 0.0, "straight kappa 0");
    }

    /* Sphere R=3: K = 1/9; plane: K = H = 0. */
    {
        double K = 0, H = 0;
        ml_status_t st = ml_surf_curvatures (spx, spy, spz, 1.0, 0.5, &K, &H);
        ASSERT_TRUE (&ctx, st == ML_SUCCESS, "sphere status");
        ASSERT_NEAR (&ctx, K, 0.1111111111111111, 1e-6, "sphere K=1/9");
        ASSERT_NEAR (&ctx, ml_fabs (H), 0.3333333333333333, 1e-6, "sphere |H|=1/3");
        st = ml_surf_curvatures (plx, ply, plz, 0.2, 0.3, &K, &H);
        ASSERT_TRUE (&ctx, st == ML_SUCCESS, "plane status");
        ASSERT_NEAR (&ctx, K, 0.0, 1e-9, "plane K=0");
        ASSERT_NEAR (&ctx, H, 0.0, 1e-9, "plane H=0");
        ASSERT_TRUE (&ctx, ml_surf_curvatures (NULL, ply, plz, 0.2, 0.3, &K, &H) == ML_ERR_INVALID_ARG, "surf null");
    }

    return ml_test_summary (&ctx);
}
