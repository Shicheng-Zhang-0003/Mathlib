/* v12R2 extended edge pins for the eleven depth-pack modules.
 * Ground truth: hand-computed / closed-form values. */
#include "test_harness.h"
#include "ml_core.h"
#include "ml_tensor.h"
#include "ml_math_extra.h"
#include "ml_bigint.h"
#include "ml_special.h"
#include "ml_linalg_extra.h"
#include "ml_calculus_extra.h"
#include "ml_stats_extra.h"
#include "ml_statistics.h"
#include "ml_graph_extra.h"
#include "ml_optim_extra.h"
#include "ml_geom_extra.h"
#include "ml_complexan.h"
#include "ml_ode_extra.h"
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
static double ef_sin (double x) { return sin (x); }
static double ef_sq (const double *x, int n, void *c) {
    (void) n;
    (void) c;
    return x [0] * x [0] + 2.0 * x [1] * x [1];
}
static void ef_ode (double t, const double *y, double *dy, int n, void *c) {
    (void) t;
    (void) c;
    for (int i = 0; i < n; i++) dy [i] = -y [i];
}
int main (void) {
    ml_test_ctx_t ctx;
    ml_test_init (&ctx, "Edge Extended v12R2");

    /* Wave 1: libm parity. */
    {
        int q = 0;
        double a = ml_remquo (17.0, 5.0, &q);
        ASSERT_NEAR (&ctx, a, 2.0, 1e-15, "remquo value");
        ASSERT_TRUE (&ctx, q == 3, "remquo quotient");
        ASSERT_TRUE (&ctx, ml_ilogb (8.0) == 3, "ilogb");
        ASSERT_NEAR (&ctx, ml_fma (2.0, 3.0, 4.0), 10.0, 0.0, "fma");
    }

    /* Wave 2: bigint 2^10 mod 1000 = 24. */
    {
        ml_bigint_t b, e, m, r;
        ml_bigint_init (&b);
        ml_bigint_init (&e);
        ml_bigint_init (&m);
        ml_bigint_init (&r);
        ml_bigint_set_u64 (&b, 2);
        ml_bigint_set_u64 (&e, 10);
        ml_bigint_set_u64 (&m, 1000);
        ASSERT_TRUE (&ctx, ml_bigint_powmod (&r, &b, &e, &m) == ML_SUCCESS, "bigint powmod st");
        ASSERT_NEAR (&ctx, ml_bigint_to_double (&r), 24.0, 0.0, "bigint powmod 2^10%1000");
    }

    /* Wave 3: special functions. */
    ASSERT_NEAR (&ctx, ml_bessel_jv (0.0, 1.0), 0.7651976865579666, 1e-12, "J0(1)");
    ASSERT_NEAR (&ctx, ml_bessel_kv (0.0, 1.0), 0.4210244382407083, 1e-12, "K0(1)");

    /* Wave 4: SVD singular values of diag(3,2). */
    {
        double A [4] = {3, 0, 0, 2}, U [4], Vt [4], S [2];
        unsigned char buf [65536];
        ml_workspace_t ws;
        ml_workspace_init (&ws, buf, sizeof (buf));
        ASSERT_TRUE (&ctx, ml_svd (ml_tensor_view (A, 2, 2), ml_tensor_view (U, 2, 2), S, ml_tensor_view (Vt, 2, 2), &ws) == ML_SUCCESS,
                     "svd st");
        ASSERT_NEAR (&ctx, S [0], 3.0, 1e-12, "svd s0");
        ASSERT_NEAR (&ctx, S [1], 2.0, 1e-12, "svd s1");
    }

    /* Wave 5: gradient + Richardson derivative. */
    {
        double x [2] = {1.0, 1.0}, g [2], err;
        ml_gradient_rn (ef_sq, x, 2, g, 1e-5, NULL);
        ASSERT_NEAR (&ctx, g [0], 2.0, 1e-6, "grad x0");
        ASSERT_NEAR (&ctx, g [1], 4.0, 1e-6, "grad x1");
        ASSERT_NEAR (&ctx, ml_derivative_richardson (ef_sin, 1.0, 1e-3, &err),
                     cos (1.0), 1e-10, "richardson dsin");
    }

    /* Wave 6: statistics. */
    {
        double d [5] = {1, 2, 3, 4, 5}, tmp [5], t, p;
        ASSERT_NEAR (&ctx, ml_quantile (d, 5, 0.5, tmp), 3.0, 0.0, "quantile");
        ASSERT_TRUE (&ctx, ml_ttest_1samp (d, 5, 3.0, &t, &p) == ML_SUCCESS, "ttest st");
        ASSERT_NEAR (&ctx, t, 0.0, 1e-12, "ttest t");
        ASSERT_NEAR (&ctx, p, 1.0, 1e-9, "ttest p");
    }

    /* Wave 7: graph. */
    {
        double W [9], D [9];
        int nx [9];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                W [i * 3 + j] = (i == j) ? 0.0 : ML_GRAPH_INF;
        W [0 * 3 + 1] = 2.0;
        W [1 * 3 + 2] = 3.0;
        ml_floyd_warshall (3, W, D, nx);
        ASSERT_NEAR (&ctx, D [0 * 3 + 2], 5.0, 1e-12, "floyd");
    }

    /* Wave 8: simplex max 3x+2y s.t. x+y<=4, x+3y<=6 -> 12. */
    {
        double c [2] = {3, 2}, A [4] = {1, 1, 1, 3}, b [2] = {4, 6}, x [2], o;
        ASSERT_TRUE (&ctx, ml_simplex_lp (c, A, b, 2, 2, x, &o) == ML_SUCCESS, "lp st");
        ASSERT_NEAR (&ctx, o, 12.0, 1e-9, "lp obj");
    }

    /* Wave 9: enclosing circle radius 2, Delaunay of square = 2 triangles. */
    {
        double x [4] = {-1, 3, 1, 1}, y [4] = {1, 1, -1, 3}, cx, cy, r;
        ml_min_enclosing_circle (x, y, 4, &cx, &cy, &r);
        ASSERT_NEAR (&ctx, r, 2.0, 1e-9, "mec r");
        {
            double sx [4] = {0, 1, 1, 0}, sy [4] = {0, 0, 1, 1};
            int tris [12];
            ASSERT_TRUE (&ctx, ml_delaunay_triangulate (sx, sy, 4, tris, 4) == 2,
                         "delaunay count");
        }
    }

    /* Wave 10: complex roots of z^2 - 3z + 2 -> {1,2}. */
    {
        cplx co [3] = {{1, 0}, {-3, 0}, {2, 0}}, roots [2];
        ASSERT_TRUE (&ctx, ml_cplx_poly_roots (co, 2, roots) == ML_SUCCESS, "cplx roots st");
        ASSERT_NEAR (&ctx, roots [0].real + roots [1].real, 3.0, 1e-9, "cplx root sum");
        ASSERT_NEAR (&ctx,
                     roots [0].real * roots [1].real - roots [0].imag * roots [1].imag,
                     2.0, 1e-9, "cplx root prod");
    }

    /* Wave 11: RK4 on y'=-y, y(0)=1 -> e^-1. */
    {
        double y0 [1] = {1.0}, out [1];
        ml_ode_rk4_sys (ef_ode, 0.0, y0, 1, 0.01, 100, NULL, out);
        ASSERT_NEAR (&ctx, out [0], exp (-1.0), 1e-9, "rk4 exp(-1)");
    }

    return ml_test_summary (&ctx);
}
