/* bench_vs_libm.c - head-to-head timing of mathlib ml_* against system libm.
 *
 * Method: for each function a deterministic input array is generated once;
 * both implementations are timed on the SAME array with direct (non-indirect)
 * calls so no call-through-pointer overhead pollutes either side.  Wall time
 * is taken with clock_gettime(CLOCK_MONOTONIC); the reported ns/op is the
 * minimum over REPS passes (warm cache, noise-robust).  A volatile sink keeps
 * the loops from being optimised away.
 *
 * Build (strict flags, same TU set as the library):
 *   gcc -std=c99 -O3 -fno-fast-math -ffp-contract=off -Iinclude/mathlib -Isrc \
 *       -o build/bench_vs_libm benchmarks/bench_vs_libm.c $SRC -lm
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdint.h>

#include "ml_core.h"
#include "ml_trig.h"
#include "ml_exp_log.h"
#include "ml_integral.h"
#include "ml_special.h"

#define NELEM 8192
#define REPS 200

static uint64_t rng_state = 0x123456789abcdefULL;
static double next_unit (void) { /* uniform [0,1) */
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return (double) (rng_state >> 11) * (1.0 / 9007199254740992.0);
}
static double next_range (double lo, double hi) { return lo + (hi - lo) * next_unit (); }
static volatile double sink = 0.0;
static double now_ns (void) {
    struct timespec ts;
    clock_gettime (CLOCK_MONOTONIC, &ts);
    return (double) ts.tv_sec * 1e9 + (double) ts.tv_nsec;
}
/* collect ratios for the summary (geomean / best / worst / slowdown buckets) */
#define MAXCASE 128
static double g_ratio [MAXCASE];
static const char *g_name [MAXCASE];
static int g_n = 0;
static void report (const char *label, double ref, double mine) {
    double r = mine / ref;
    if (g_n < MAXCASE) {
        g_ratio [g_n] = r;
        g_name [g_n] = label;
        g_n++;
    }
    printf ("%-9s %10.2f %10.2f %8.2f\n", label, ref, mine, r);
}
/* ---- unary: libm is the reference (ref), mathlib is mine ---- */
#define DEF_U(mine, ref, LO, HI, label)                              \
    static void bench_##mine (void) {                                \
        static double x [NELEM];                                     \
        int i, r, k;                                                 \
        double acc, t0, t1, bm = 1e300, br = 1e300;                  \
        for (i = 0; i < NELEM; i++) x [i] = next_range ((LO), (HI)); \
        for (r = 0; r < REPS; r++) {                                 \
            acc = 0.0;                                               \
            t0 = now_ns ();                                          \
            for (k = 0; k < NELEM; k++) acc += mine (x [k]);         \
            t1 = now_ns ();                                          \
            sink += acc;                                             \
            if (t1 - t0 < bm) bm = t1 - t0;                          \
            acc = 0.0;                                               \
            t0 = now_ns ();                                          \
            for (k = 0; k < NELEM; k++) acc += ref (x [k]);          \
            t1 = now_ns ();                                          \
            sink += acc;                                             \
            if (t1 - t0 < br) br = t1 - t0;                          \
        }                                                            \
        report (label, br / NELEM, bm / NELEM);                      \
    }

/* ---- binary ---- */
#define DEF_B(mine, ref, ALO, AHI, BLO, BHI, label)                 \
    static void bench_##mine (void) {                               \
        static double a [NELEM], b [NELEM];                         \
        int i, r, k;                                                \
        double acc, t0, t1, bm = 1e300, br = 1e300;                 \
        for (i = 0; i < NELEM; i++) {                               \
            a [i] = next_range ((ALO), (AHI));                      \
            b [i] = next_range ((BLO), (BHI));                      \
        }                                                           \
        for (r = 0; r < REPS; r++) {                                \
            acc = 0.0;                                              \
            t0 = now_ns ();                                         \
            for (k = 0; k < NELEM; k++) acc += mine (a [k], b [k]); \
            t1 = now_ns ();                                         \
            sink += acc;                                            \
            if (t1 - t0 < bm) bm = t1 - t0;                         \
            acc = 0.0;                                              \
            t0 = now_ns ();                                         \
            for (k = 0; k < NELEM; k++) acc += ref (a [k], b [k]);  \
            t1 = now_ns ();                                         \
            sink += acc;                                            \
            if (t1 - t0 < br) br = t1 - t0;                         \
        }                                                           \
        report (label, br / NELEM, bm / NELEM);                     \
    }

DEF_U (ml_sqrt, sqrt, 1.0, 1.0e6, "sqrt")
DEF_U (ml_cbrt, cbrt, 1.0e-3, 1.0e6, "cbrt")
DEF_U (ml_exp, exp, -700.0, 700.0, "exp")
DEF_U (ml_exp2, exp2, -1000.0, 1000.0, "exp2")
DEF_U (ml_exp10, exp10, -300.0, 300.0, "exp10")
DEF_U (ml_expm1, expm1, -30.0, 30.0, "expm1")
DEF_U (ml_log, log, 1.0e-300, 1.0e300, "log")
DEF_U (ml_log2, log2, 1.0e-300, 1.0e300, "log2")
DEF_U (ml_log10, log10, 1.0e-300, 1.0e300, "log10")
DEF_U (ml_log1p, log1p, -0.9, 1.0e6, "log1p")
DEF_U (ml_sin, sin, -100.0, 100.0, "sin")
DEF_U (ml_cos, cos, -100.0, 100.0, "cos")
DEF_U (ml_tan, tan, -1.5, 1.5, "tan")
DEF_U (ml_asin, asin, -0.999, 0.999, "asin")
DEF_U (ml_acos, acos, -0.999, 0.999, "acos")
DEF_U (ml_atan, atan, -1000.0, 1000.0, "atan")
DEF_U (ml_sinh, sinh, -30.0, 30.0, "sinh")
DEF_U (ml_cosh, cosh, -30.0, 30.0, "cosh")
DEF_U (ml_tanh, tanh, -30.0, 30.0, "tanh")
DEF_U (ml_asinh, asinh, -100.0, 100.0, "asinh")
DEF_U (ml_acosh, acosh, 1.0, 100.0, "acosh")
DEF_U (ml_atanh, atanh, -0.999, 0.999, "atanh")
DEF_U (ml_erf, erf, -6.0, 6.0, "erf")
DEF_U (ml_erfc, erfc, -6.0, 6.0, "erfc")
DEF_U (ml_gamma_new, tgamma, 0.1, 20.0, "tgamma")
DEF_U (ml_lgamma, lgamma, 0.1, 20.0, "lgamma")
DEF_U (ml_bessel_j0, j0, 0.1, 50.0, "j0")
DEF_U (ml_bessel_j1, j1, 0.1, 50.0, "j1")
DEF_U (ml_bessel_y0, y0, 0.1, 50.0, "y0")
DEF_U (ml_bessel_y1, y1, 0.1, 50.0, "y1")

DEF_B (ml_pow, pow, 0.1, 10.0, -5.0, 5.0, "pow")
DEF_B (ml_atan2, atan2, -100.0, 100.0, -100.0, 100.0, "atan2")
int main (void) {
    printf ("mathlib vs system libm (glibc)  --  ns/op, min of %d passes over %d elems\n",
            REPS, NELEM);
    printf ("=================================================================\n");
    printf ("%-9s %10s %10s %8s\n", "func", "libm", "mathlib", "ratio");
    printf ("-----------------------------------------------------------------\n");
    bench_ml_sqrt ();
    bench_ml_cbrt ();
    bench_ml_exp ();
    bench_ml_exp2 ();
    bench_ml_exp10 ();
    bench_ml_expm1 ();
    bench_ml_log ();
    bench_ml_log2 ();
    bench_ml_log10 ();
    bench_ml_log1p ();
    bench_ml_sin ();
    bench_ml_cos ();
    bench_ml_tan ();
    bench_ml_asin ();
    bench_ml_acos ();
    bench_ml_atan ();
    bench_ml_sinh ();
    bench_ml_cosh ();
    bench_ml_tanh ();
    bench_ml_asinh ();
    bench_ml_acosh ();
    bench_ml_atanh ();
    bench_ml_erf ();
    bench_ml_erfc ();
    bench_ml_gamma_new ();
    bench_ml_lgamma ();
    bench_ml_bessel_j0 ();
    bench_ml_bessel_j1 ();
    bench_ml_bessel_y0 ();
    bench_ml_bessel_y1 ();
    bench_ml_pow ();
    bench_ml_atan2 ();
    printf ("=================================================================\n");

    /* summary: geomean/median slowdown, best & worst.  sqrt is excluded
     * because glibc's sqrt is a single hardware instruction the compiler
     * inlines, an apples-to-oranges baseline; see -fno-builtin control. */
    {
        double sumln = 0.0, best = 1e300, worst = 0.0;
        int n = 0, le2 = 0, gt10 = 0, i;
        const char *bestn = "", *worstn = "";
        for (i = 0; i < g_n; i++) {
            if (g_ratio [i] <= 0.0) continue;
            if (!(g_name [i][0] == 's' && g_name [i][1] == 'q')) {
                sumln += log (g_ratio [i]);
                n++;
            }
            if (g_ratio [i] < best) {
                best = g_ratio [i];
                bestn = g_name [i];
            }
            if (g_ratio [i] > worst) {
                worst = g_ratio [i];
                worstn = g_name [i];
            }
            if (g_ratio [i] <= 2.0) le2++;
            if (g_ratio [i] > 10.0) gt10++;
        }
        printf ("geomean mathlib/libm slowdown (excl. sqrt): %.2fx over %d funcs\n",
                exp (sumln / (double) n), n);
        printf ("best %s %.2fx | worst %s %.2fx | <=2x: %d | >10x: %d\n",
                bestn, best, worstn, worst, le2, gt10);
    }
    printf ("ratio = mathlib/libm  (<1 faster, >1 slower); sink=%g\n", sink);
    return 0;
}
