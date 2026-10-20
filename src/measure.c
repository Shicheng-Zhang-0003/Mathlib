#include "ml_compiler.h"
#include "ml_measure.h"
#include "ml_core.h"
#include <stddef.h>
#include <stdlib.h>
/* Measure-theoretic instance verifiers on uniform samples. Long-double
 * accumulation; verdicts fail closed. */

ML_API double ml_lebesgue_dirichlet (double a, double b) {
    /* 1_Q = 0 a.e.: the integral over any valid interval is exactly 0. */
    if (ml_isnan (a) || ml_isnan (b)) return ml_make_nan ();
    if (!ml_isfinite (a) || !ml_isfinite (b)) return ml_make_nan ();
    if (a > b) return ml_make_nan ();
    return 0.0;
}
ML_API double ml_simple_integral (const double *values, const double *measures, int n) {
    long double s = 0.0L;
    int i;
    double r;
    if (!values || !measures || n <= 0 || n > 1000000) return ml_make_nan ();
    for (i = 0; i < n; i++) {
        if (ml_isnan (values [i]) || ml_isnan (measures [i])) return ml_make_nan ();
        if (!ml_isfinite (values [i]) || !ml_isfinite (measures [i])) return ml_make_nan ();
        if (measures [i] < 0.0) return ml_make_nan ();
        s += (long double) values [i] * (long double) measures [i];
    }
    r = (double) s;
    return ml_isfinite (r) ? r : ml_make_nan ();
}
ML_API double ml_lp_norm_data (const double *f, int n, double dx, double p) {
    /* (sum |f|^p dx)^{1/p}, p >= 1. */
    long double s = 0.0L;
    int i;
    double r;
    if (!f || n <= 0 || n > 1000000) return ml_make_nan ();
    if (ml_isnan (dx) || ml_isnan (p)) return ml_make_nan ();
    if (!ml_isfinite (dx) || !(dx > 0.0)) return ml_make_nan ();
    if (!ml_isfinite (p) || !(p >= 1.0)) return ml_make_nan ();
    for (i = 0; i < n; i++) {
        long double v;
        if (ml_isnan (f [i]) || !ml_isfinite (f [i])) return ml_make_nan ();
        v = (long double) f [i];
        if (v < 0.0L) v = -v;
        s += __builtin_powl (v, p);
    }
    s *= (long double) dx;
    r = (double) __builtin_powl (s, 1.0L / p);
    return ml_isfinite (r) ? r : ml_make_nan ();
}
ML_API double ml_seq_lp_norm (const double *a, int n, double p) {
    return ml_lp_norm_data (a, n, 1.0, p);
}
ML_API int ml_holder_check (const double *f, const double *g, int n,
                            double dx, double p, double tol) {
    /* sum|fg|dx <= ||f||_p ||g||_q, 1/p + 1/q = 1. */
    long double lhs = 0.0L;
    double nf;
    double ng;
    double q;
    int i;
    if (!f || !g || n <= 0 || n > 1000000) return 0;
    if (ml_isnan (dx) || !ml_isfinite (dx) || !(dx > 0.0)) return 0;
    if (ml_isnan (p) || !ml_isfinite (p) || !(p > 1.0)) return 0;
    if (!(tol > 0.0) || !ml_isfinite (tol)) tol = 1e-9;
    q = p / (p - 1.0);
    for (i = 0; i < n; i++) {
        long double a;
        long double b;
        if (ml_isnan (f [i]) || ml_isnan (g [i])) return 0;
        if (!ml_isfinite (f [i]) || !ml_isfinite (g [i])) return 0;
        a = (long double) f [i];
        b = (long double) g [i];
        if (a < 0.0L) a = -a;
        if (b < 0.0L) b = -b;
        lhs += a * b;
    }
    lhs *= (long double) dx;
    nf = ml_lp_norm_data (f, n, dx, p);
    ng = ml_lp_norm_data (g, n, dx, q);
    if (!ml_isfinite (nf) || !ml_isfinite (ng)) return 0;
    {
        long double rhs = (long double) nf * (long double) ng;
        long double sl = tol * (1.0L + rhs);
        return (lhs <= rhs + sl) ? 1 : 0;
    }
}
ML_API int ml_minkowski_check (const double *f, const double *g, int n,
                               double dx, double p, double tol) {
    /* ||f+g||_p <= ||f||_p + ||g||_p, p >= 1. */
    double *h = NULL;
    double lhs;
    double nf;
    double ng;
    int i;
    int ok;
    /* Stack buffer covers n <= 1024 (8KB); larger goes heap-per-call. */
    double stack [1024];
    if (!f || !g || n <= 0 || n > 1000000) return 0;
    if (ml_isnan (dx) || !ml_isfinite (dx) || !(dx > 0.0)) return 0;
    if (ml_isnan (p) || !ml_isfinite (p) || !(p >= 1.0)) return 0;
    if (!(tol > 0.0) || !ml_isfinite (tol)) tol = 1e-9;
    if (n <= 1024) {
        h = stack;
    } else {
        if ((size_t) n > (size_t) -1 / sizeof (double)) return 0;
        h = (double *) malloc ((size_t) n * sizeof (double));
        if (!h) return 0;
    }
    for (i = 0; i < n; i++) {
        if (ml_isnan (f [i]) || ml_isnan (g [i])) {
            if (h != stack) free (h);
            return 0;
        }
        if (!ml_isfinite (f [i]) || !ml_isfinite (g [i])) {
            if (h != stack) free (h);
            return 0;
        }
        h [i] = f [i] + g [i];
    }
    lhs = ml_lp_norm_data (h, n, dx, p);
    nf = ml_lp_norm_data (f, n, dx, p);
    ng = ml_lp_norm_data (g, n, dx, p);
    if (h != stack) free (h);
    if (!ml_isfinite (lhs) || !ml_isfinite (nf) || !ml_isfinite (ng)) return 0;
    ok = (lhs <= nf + ng + tol * (1.0 + nf + ng)) ? 1 : 0;
    return ok;
}
ML_API double ml_markov_bound (double mean_abs, double t) {
    /* P(|X| >= t) <= E|X|/t. */
    if (ml_isnan (mean_abs) || ml_isnan (t)) return ml_make_nan ();
    if (!ml_isfinite (mean_abs) || !ml_isfinite (t)) return ml_make_nan ();
    if (!(mean_abs >= 0.0) || !(t > 0.0)) return ml_make_nan ();
    return mean_abs / t;
}
