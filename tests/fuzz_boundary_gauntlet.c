#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include <time.h>
#include "ml_core.h"
#include "bitwise_fp.h"
#include "ml_trig.h"
#include "ml_exp_log.h"
#include "ml_fixed_point.h"
#include "ml_tensor.h"
#include "ml_linalg.h"
#include "ml_types.h"
#include "fft.h"
#include "ml_complex.h"
#include "fast_math.h"
#include "ml_contest.h"
#include "ml_numbertheory.h"
#include "ml_analytic_nt.h"
#include "ml_geometry.h"
#include "ml_graph.h"
#include "ml_modeling.h"
#include "ml_topology.h"
#include "ml_diffgeom.h"
#include "ml_algebraic.h"
#include "ml_funcanal.h"

int passed = 0;
int failed = 0;

#define CHECK(cond, msg)                                    \
    do {                                                    \
        if (cond) {                                         \
            passed++;                                       \
        } else {                                            \
            failed++;                                       \
            printf ("FAIL: %s (Line %d)\n", msg, __LINE__); \
        }                                                   \
    } while (0)

#define CHECK_NEAR(a, b, eps, msg) CHECK (ml_fabs ((double) (a) - (double) (b)) < (eps), msg)
void test_trig_boundaries () {
    printf ("--- Trig Exact Boundaries ---\n");
    CHECK_NEAR (ml_sin (0.0), 0.0, 1e-15, "sin(0)");
    CHECK_NEAR (ml_sin (ML_PI / 2.0), 1.0, 1e-15, "sin(pi/2)");
    CHECK_NEAR (ml_sin (ML_PI), 0.0, 1e-14, "sin(pi)");
    CHECK_NEAR (ml_sin (3.0 * ML_PI / 2.0), -1.0, 1e-14, "sin(3pi/2)");
    CHECK_NEAR (ml_sin (2.0 * ML_PI), 0.0, 1e-14, "sin(2pi)");
    CHECK_NEAR (ml_cos (0.0), 1.0, 1e-15, "cos(0)");
    CHECK_NEAR (ml_cos (ML_PI / 2.0), 0.0, 1e-15, "cos(pi/2)");
    CHECK_NEAR (ml_cos (ML_PI), -1.0, 1e-14, "cos(pi)");

    // Boundary crossing: 1 ULP away from pi/2
    double eps = 2.220446049250313e-16;
    CHECK_NEAR (ml_sin (ML_PI / 2.0 - eps), ml_cos (eps), 1e-12, "sin(pi/2 - eps) == cos(eps)");
    CHECK_NEAR (ml_sin (ML_PI / 2.0 + eps), ml_cos (eps), 1e-12, "sin(pi/2 + eps) == cos(eps)");
    CHECK_NEAR (ml_sin (-ML_PI / 2.0 - eps), -ml_cos (eps), 1e-12, "sin(-pi/2 - eps) == -cos(eps)");
}
void test_exp_log_boundaries () {
    printf ("--- Exp/Log Exact Boundaries ---\n");
    CHECK_NEAR (ml_log (1.0), 0.0, 1e-15, "log(1)");
    CHECK_NEAR (ml_log (ML_E), 1.0, 1e-14, "log(e)");
    CHECK_NEAR (ml_exp (0.0), 1.0, 1e-15, "exp(0)");
    CHECK_NEAR (ml_exp (1.0), ML_E, 1e-14, "exp(1)");
    CHECK (ml_isinf (ml_log (0.0)) && ml_log (0.0) < 0, "log(0) == -inf");
    CHECK (ml_isinf (ml_exp (1000.0)) && ml_exp (1000.0) > 0, "exp(1000) == +inf");
    CHECK (ml_exp (-1000.0) < 5e-324, "exp(-1000) underflows to zero or subnormal");
}
void test_fixed_point_saturation () {
    printf ("--- Fixed-Point Saturation ---\n");
    ml_q16_16_t max_val = 2147483647;
    ml_q16_16_t min_val = (-2147483647 - 1);
    CHECK (ml_fixed_mul (max_val, max_val) == 2147483647, "MAX * MAX saturates to MAX");
    CHECK (ml_fixed_mul (min_val, max_val) == (-2147483647 - 1), "MIN * MAX saturates to MIN");
    CHECK (ml_fixed_mul (min_val, min_val) == 2147483647, "MIN * MIN saturates to MAX");
    ml_q16_16_t two = 2 << 16;
    ml_q16_16_t three = 3 << 16;
    CHECK (ml_fixed_mul (two, three) == (6 << 16), "2.0 * 3.0 == 6.0");
}
void test_tensor_hilbert () {
    printf ("--- Tensor Ill-Conditioned (Hilbert) ---\n");
    char scratchpad [4096]; // FIXED: 1MB stack bomb reduced to 4KB
    ml_workspace_t ws;
    ml_workspace_init (&ws, scratchpad, sizeof (scratchpad));
    int n = 4;
    double A_data [16];
    double b_data [4] = {1.0, 1.0, 1.0, 1.0};
    double x_data [4] = {0};
    double residual [4] = {0};

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A_data [i * n + j] = 1.0 / (i + j + 1.0);
        }
    }
    ml_tensor_view_t A_view = ml_tensor_view (A_data, n, n);
    int status = ml_solve (A_view, b_data, x_data, &ws);
    CHECK (status == ML_SUCCESS, "Hilbert solve status");

    // FIXED: Actually validate the solution using residual ||Ax - b||
    ml_matvec (A_view, x_data, residual);
    double max_err = 0.0;
    for (int i = 0; i < n; i++) {
        double err = ml_fabs (residual [i] - b_data [i]);
        if (err > max_err) max_err = err;
    }
    // Hilbert n=4 condition number is ~15,000. 1e-10 is a strict but fair bound.
    CHECK (max_err < 1e-10, "Hilbert residual ||Ax - b|| < 1e-10");
}
void test_fft_parseval () {
    printf ("--- FFT Parseval's Theorem (Energy Conservation) ---\n");
    srand (42); // FIXED: Deterministic seed for exact reproduction
    int n = 64;
    cplx time_domain [64];
    cplx freq_domain [64];
    double time_energy = 0.0;
    for (int i = 0; i < n; i++) {
        double val = (double) rand () / RAND_MAX * 10.0 - 5.0;
        time_domain [i] = (cplx) {val, 0.0};
        time_energy += val * val;
    }
    for (int i = 0; i < n; i++) freq_domain [i] = time_domain [i];
    ml_fft_execute (freq_domain, n);
    double freq_energy = 0.0;
    for (int i = 0; i < n; i++) {
        freq_energy += (freq_domain [i].real * freq_domain [i].real + freq_domain [i].imag * freq_domain [i].imag);
    }
    freq_energy /= n;
    // FIXED: Tightened from 1e-6 (useless) to 1e-14 (physical limit of double precision)
    CHECK (ml_fabs (time_energy - freq_energy) < time_energy * 1e-13, "Parseval's Theorem (Time Energy == Freq Energy)");
}
void test_contest_boundaries () {
    printf ("--- Contest Exact Boundaries ---\n");
    CHECK (ml_narayana (5, 1) + ml_narayana (5, 2) + ml_narayana (5, 3) + ml_narayana (5, 4) + ml_narayana (5, 5) == 42, "narayana row 5 == catalan 5");
    CHECK (ml_hook_dimension ((int []) {3, 2, 1}, 3) == 16, "hook 321 == 16");
    CHECK (ml_sn_class_size ((int []) {2, 2}, 2) == 3, "class 2+2 == 3");
    CHECK (ml_lte_vp (10, 4, 6, 3) == 2, "LTE textbook v3(10^6-4^6)");
    CHECK (ml_floor_sum_u64 (10, 1, 3, 2) == 155, "floor_sum closed form");
    CHECK (ml_sqrt_cf_period (23) == 4, "cf period sqrt23");
    {
        uint64_t x = 0, y = 0;
        CHECK (ml_pell_fund (2, &x, &y) == ML_SUCCESS && x == 3 && y == 2, "pell D=2");
    }
    CHECK (ml_collatz_steps (27) == 111, "collatz 27");
}
void test_geometry_boundaries () {
    printf ("--- Geometry Exact Boundaries ---\n");
    CHECK_NEAR (ml_tri_area_coords (0, 0, 4, 0, 0, 3), 6.0, 1e-12, "3-4-5 area");
    CHECK_NEAR (ml_tri_circumradius (3, 4, 5), 2.5, 1e-12, "R == abc/4K");
    CHECK_NEAR (ml_tri_inradius (3, 4, 5), 1.0, 1e-12, "r == K/s");
    CHECK (ml_euler_line_check (0, 0, 4, 0, 0, 3, 1e-9) == 1, "euler line holds");
    CHECK_NEAR (ml_power_of_point (1.0, 2.0, 1.0, 2.0, 3.0), -9.0, 1e-12, "power at center");
    CHECK_NEAR (ml_clock_angle (3, 0), 90.0, 1e-12, "clock 3:00");
    {
        double xs [] = {0, 1, 1, 0}, ys [] = {0, 0, 1, 1};
        int out [4] = {0, 0, 0, 0};
        CHECK (ml_convex_hull (xs, ys, 4, out, 4) == 4, "hull square");
    }
}
void test_graph_boundaries () {
    printf ("--- Graph Exact Boundaries ---\n");
    {
        int eu [] = {0, 1}, ev [] = {1, 2};
        double w [] = {1.0, 2.0};
        double d [3] = {0, 0, 0};
        double tot = 0;
        CHECK (ml_dijkstra (3, eu, ev, w, 2, 0, d) == ML_SUCCESS, "dijkstra status");
        CHECK_NEAR (d [2], 3.0, 1e-12, "dijkstra path weight");
        CHECK (ml_kruskal_mst (3, eu, ev, w, 2, &tot) == ML_SUCCESS && tot == 3.0, "kruskal path");
    }
    {
        int eu [] = {0, 0, 1, 2}, ev [] = {1, 2, 3, 3};
        int out [4] = {0, 0, 0, 0};
        int c = ml_topo_sort (4, eu, ev, 4, out);
        CHECK (c == 4 && out [0] == 0 && out [3] == 3, "topo diamond");
    }
    {
        int e3 [] = {0, 1, 2}, f3 [] = {1, 2, 0};
        CHECK (ml_is_bipartite (3, e3, f3, 3) == 0, "odd cycle not bipartite");
    }
    {
        int p [4], r [4];
        ml_dsu_make (p, r, 4);
        ml_dsu_union (p, r, 4, 0, 3);
        CHECK (ml_dsu_find (p, 4, 0) == ml_dsu_find (p, 4, 3), "dsu union-find");
    }
}
static double bnd_banach_t (double x) {
    return 0.5 * x + 0.25;
}
void test_modeling_boundaries () {
    printf ("--- Modeling Closed Forms ---\n");
    {
        double P [] = {0.7, 0.3, 0.4, 0.6};
        double pi [2] = {0, 0};
        CHECK (ml_markov_stationary (P, pi, 2, 10000, 1e-12) == ML_SUCCESS, "markov status");
        CHECK_NEAR (pi [0], 0.5714285714285714, 1e-9, "stationary 4/7");
    }
    CHECK_NEAR (ml_sir_final_size (3.0), 0.9404797907, 1e-7, "SIR Z(3)");
    CHECK_NEAR (ml_mm1_wait (2, 3, 1.0), 0.36787944117144233, 1e-12, "M/M/1 tail e^-1");
    {
        double a1 = 0, b1 = 0;
        CHECK (ml_beta_binomial_post (2, 2, 7, 10, &a1, &b1) == ML_SUCCESS && a1 == 9.0 && b1 == 5.0, "conjugate exact");
    }
    {
        double xs = 0;
        CHECK (ml_banach_fixed_point (bnd_banach_t, 0.0, 0.5, 1e-12, 10000, &xs, NULL) == ML_SUCCESS && ml_fabs (xs - 0.5) < 1e-9, "banach 1/2");
    }
}
void test_doctorate_boundaries () {
    printf ("--- Doctorate Exact Boundaries ---\n");
    {
        int v [] = {0, 1, 2, 3, 0, 1, 0, 2, 0, 3, 1, 2, 1, 3, 2, 3,
                    0, 1, 2, 0, 1, 3, 0, 2, 3, 1, 2, 3};
        int off [] = {0, 1, 2, 3, 4, 6, 8, 10, 12, 14, 16, 19, 22, 25, 28};
        int dim [] = {0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2};
        int b [4] = {0, 0, 0, 0};
        CHECK (ml_simplicial_betti_f2 (4, v, off, dim, 14, b, 3) == ML_SUCCESS, "tet betti status");
        CHECK (b [0] == 1 && b [1] == 0 && b [2] == 1, "tet betti (1,0,1)");
    }
    {
        long long d = 0;
        CHECK (ml_bareiss_det ((long long []) {1, 2, 3, 4}, 2, &d) == ML_SUCCESS && d == -2, "bareiss det");
        CHECK (ml_sturm_count ((long long []) {-6, 11, -6, 1}, 3, 0.0, 4.0) == 3, "sturm cubic");
    }
    {
        double R [9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double z [3] = {0, 0, 1};
        CHECK (ml_rodrigues_exp (z, 0.0, R) == ML_SUCCESS, "rodrigues status");
        CHECK (R [0] == 1.0 && R [4] == 1.0 && R [8] == 1.0 && R [1] == 0.0, "rodrigues identity");
    }
    CHECK_NEAR (ml_ei (1.0), 1.8951178163559368, 1e-9, "Ei(1)");
    CHECK_NEAR (ml_li (1000000.0), 78627.5492, 0.01, "li(1e6)");
    {
        CHECK (ml_fft_is_supported (1024) == 1, "fft pow2 supported");
        CHECK (ml_fft_is_supported (1000) == 0, "fft non-pow2 refused");
        CHECK (ml_fft_is_supported (1 << 25) == 0, "fft 2^25 refused (DoS cap)");
    }
}
int main () {
    printf ("=========================================================\n");
    printf ("   MATHLIB v12R2: BOUNDARY & INVARIANT GAUNTLET\n");
    printf ("=========================================================\n");
    test_trig_boundaries ();
    test_exp_log_boundaries ();
    test_fixed_point_saturation ();
    test_tensor_hilbert ();
    test_fft_parseval ();
    test_contest_boundaries ();
    test_geometry_boundaries ();
    test_graph_boundaries ();
    test_modeling_boundaries ();
    test_doctorate_boundaries ();
    printf ("\n=========================================================\n");
    printf ("BOUNDARY SUMMARY: %d passed, %d failed\n", passed, failed);
    printf ("=========================================================\n");
    return failed > 0 ? 1 : 0;
}
