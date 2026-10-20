#include "ml_compiler.h"
#include "ml_diffgeom.h"
#include "ml_core.h"
#include "ml_numerical.h"
#include <stddef.h>
/* Curves and parametric surfaces via central differences (LD). */

ML_API double ml_curvature_graph (ml_func_t f, double x) {
    /* kappa = |f''|/(1+f'^2)^{3/2}, staged derivative kernels. */
    double h;
    double fp;
    double fpp;
    long double den;
    double o;
    if (!f) return ml_make_nan ();
    if (ml_isnan (x) || !ml_isfinite (x)) return ml_make_nan ();
    h = 1.49e-8 * (1.0 + ml_fabs (x));
    fp = ml_derivative (f, x, h);
    fpp = ml_second_derivative (f, x, h);
    if (ml_isnan (fp) || ml_isnan (fpp)) return ml_make_nan ();
    if (!ml_isfinite (fp) || !ml_isfinite (fpp)) return ml_make_nan ();
    den = 1.0L + (long double) fp * (long double) fp;
    o = (double) (fabsl ((long double) fpp) / (den * __builtin_sqrtl (den)));
    return ml_isfinite (o) ? o : ml_make_nan ();
}
static long double mld_d1 (ml_func_t f, double t, long double h) {
    /* First derivative, 4th-order central: (-f2+8f1-8fm1+fm2)/12h. */
    long double f1 = (long double) f (t + (double) h);
    long double fm1 = (long double) f (t - (double) h);
    long double f2 = (long double) f (t + 2.0 * (double) h);
    long double fm2 = (long double) f (t - 2.0 * (double) h);
    return (-f2 + 8.0L * f1 - 8.0L * fm1 + fm2) / (12.0L * h);
}
static long double mld_d2 (ml_func_t f, double t, long double h) {
    /* Second derivative, 4th-order central. */
    long double f1 = (long double) f (t + (double) h);
    long double fm1 = (long double) f (t - (double) h);
    long double f2 = (long double) f (t + 2.0 * (double) h);
    long double fm2 = (long double) f (t - 2.0 * (double) h);
    long double f0 = (long double) f (t);
    return (-f2 + 16.0L * f1 - 30.0L * f0 + 16.0L * fm1 - fm2) / (12.0L * h * h);
}
static long double mld_d3 (ml_func_t f, double t, long double h) {
    /* Third derivative, 2nd-order central 4-point stencil. */
    long double f1 = (long double) f (t + (double) h);
    long double fm1 = (long double) f (t - (double) h);
    long double f2 = (long double) f (t + 2.0 * (double) h);
    long double fm2 = (long double) f (t - 2.0 * (double) h);
    return (f2 - 2.0L * f1 + 2.0L * fm1 - fm2) / (2.0L * h * h * h);
}
static long double mld_step1 (double t) {
    /* 4th-order first derivative: truncation h^4 vs rounding eps/h. */
    return 1e-5L * (1.0L + fabsl ((long double) t));
}
static long double mld_step2 (double t) {
    /* 4th-order second derivative: 64-eval coefficient weight needs
     * h >> 1e-8 in double precision (h = 2e-8 loses ~60% to rounding). */
    return 3e-4L * (1.0L + fabsl ((long double) t));
}
static long double mld_step3 (double t) {
    /* 2nd-order third derivative: truncation h^2 vs rounding eps/h^3. */
    return 1e-3L * (1.0L + fabsl ((long double) t));
}
static int mld_finite_f (ml_func_t f, double t, long double h) {
    double vs [5] = {
        f (t), f (t + (double) h), f (t - (double) h),
        f (t + 2.0 * (double) h), f (t - 2.0 * (double) h)};
    int i;
    for (i = 0; i < 5; i++) {
        if (ml_isnan (vs [i]) || !ml_isfinite (vs [i])) return 0;
    }
    return 1;
}
ML_API double ml_curvature_param (ml_func_t fx, ml_func_t fy, double t) {
    /* kappa = |x'y''-y'x''|/(x'^2+y'^2)^{3/2}. */
    long double h;
    long double x1, y1, x2, y2;
    long double sp;
    long double num;
    double o;
    if (!fx || !fy) return ml_make_nan ();
    if (ml_isnan (t) || !ml_isfinite (t)) return ml_make_nan ();
    h = mld_step2 (t);
    if (!mld_finite_f (fx, t, h) || !mld_finite_f (fy, t, h)) return ml_make_nan ();
    x1 = mld_d1 (fx, t, mld_step1 (t));
    y1 = mld_d1 (fy, t, mld_step1 (t));
    x2 = mld_d2 (fx, t, h);
    y2 = mld_d2 (fy, t, h);
    sp = x1 * x1 + y1 * y1;
    if (!(sp > 0.0L) || !ml_isfinite ((double) sp)) return ml_make_nan ();
    num = x1 * y2 - y1 * x2;
    if (num < 0.0L) num = -num;
    o = (double) (num / (sp * __builtin_sqrtl (sp)));
    return ml_isfinite (o) ? o : ml_make_nan ();
}
ML_API ml_status_t ml_frenet_frame (ml_func_t fx, ml_func_t fy, ml_func_t fz,
                                    double t, double T [3], double N [3],
                                    double B [3], double *kappa, double *tau) {
    /* T = r'/|r'|; B = r'x r''/|..|; N = B x T; kappa = |T'|/|r'|
     * via |r'x r''|/|r'|^3; tau = (r'x r'') . r'''/|r'x r''|^2. */
    long double h;
    long double r1 [3], r2 [3], r3 [3];
    long double sp;
    long double c [3];
    long double cn;
    long double kap;
    long double tor;
    long double tx [3];
    long double bx [3];
    long double nx [3];
    int i;
    if (ML_UNLIKELY (!fx || !fy || !fz || !T || !N || !B || !kappa || !tau))
        return ML_ERR_INVALID_ARG;
    if (ml_isnan (t) || !ml_isfinite (t)) {
        for (i = 0; i < 3; i++) T [i] = N [i] = B [i] = ml_make_nan ();
        *kappa = ml_make_nan ();
        *tau = ml_make_nan ();
        return ML_ERR_NAN_INPUT;
    }
    h = mld_step3 (t);
    if (!mld_finite_f (fx, t, h) || !mld_finite_f (fy, t, h) || !mld_finite_f (fz, t, h)) {
        for (i = 0; i < 3; i++) T [i] = N [i] = B [i] = ml_make_nan ();
        *kappa = ml_make_nan ();
        *tau = ml_make_nan ();
        return ML_ERR_NAN_INPUT;
    }
    r1 [0] = mld_d1 (fx, t, mld_step1 (t));
    r1 [1] = mld_d1 (fy, t, mld_step1 (t));
    r1 [2] = mld_d1 (fz, t, mld_step1 (t));
    r2 [0] = mld_d2 (fx, t, mld_step2 (t));
    r2 [1] = mld_d2 (fy, t, mld_step2 (t));
    r2 [2] = mld_d2 (fz, t, mld_step2 (t));
    r3 [0] = mld_d3 (fx, t, h);
    r3 [1] = mld_d3 (fy, t, h);
    r3 [2] = mld_d3 (fz, t, h);
    sp = r1 [0] * r1 [0] + r1 [1] * r1 [1] + r1 [2] * r1 [2];
    if (!(sp > 0.0L) || !ml_isfinite ((double) sp)) {
        for (i = 0; i < 3; i++) T [i] = N [i] = B [i] = ml_make_nan ();
        *kappa = ml_make_nan ();
        *tau = ml_make_nan ();
        return ML_ERR_SINGULAR;
    }
    c [0] = r1 [1] * r2 [2] - r1 [2] * r2 [1];
    c [1] = r1 [2] * r2 [0] - r1 [0] * r2 [2];
    c [2] = r1 [0] * r2 [1] - r1 [1] * r2 [0];
    cn = __builtin_sqrtl (c [0] * c [0] + c [1] * c [1] + c [2] * c [2]);
    if (!(cn > 0.0L) || !ml_isfinite ((double) cn)) {
        /* Straight motion: tangent exists, normal plane undefined. */
        long double n0 = __builtin_sqrtl (sp);
        for (i = 0; i < 3; i++) {
            T [i] = (double) (r1 [i] / n0);
            N [i] = ml_make_nan ();
            B [i] = ml_make_nan ();
        }
        *kappa = 0.0;
        *tau = ml_make_nan ();
        return ML_ERR_SINGULAR;
    }
    kap = cn / (sp * __builtin_sqrtl (sp));
    tor = (c [0] * r3 [0] + c [1] * r3 [1] + c [2] * r3 [2]) / (cn * cn);
    for (i = 0; i < 3; i++) {
        tx [i] = r1 [i] / __builtin_sqrtl (sp);
        bx [i] = c [i] / cn;
    }
    nx [0] = bx [1] * tx [2] - bx [2] * tx [1];
    nx [1] = bx [2] * tx [0] - bx [0] * tx [2];
    nx [2] = bx [0] * tx [1] - bx [1] * tx [0];
    for (i = 0; i < 3; i++) {
        T [i] = (double) tx [i];
        N [i] = (double) nx [i];
        B [i] = (double) bx [i];
    }
    *kappa = ml_isfinite ((double) kap) ? (double) kap : ml_make_nan ();
    *tau = ml_isfinite ((double) tor) ? (double) tor : ml_make_nan ();
    if (!ml_isfinite (*kappa)) return ML_ERR_OVERFLOW;
    return ML_SUCCESS;
}
ML_API ml_status_t ml_surf_curvatures (ml_surf_func_t rx, ml_surf_func_t ry,
                                       ml_surf_func_t rz, double u, double v,
                                       double *K, double *H) {
    /* Central differences for r_u, r_v, r_uu, r_uv, r_vv; E,F,G / L,M,N;
     * K = (LN-M^2)/(EG-F^2), H = (eG-2fF+gE)/2(EG-F^2). */
    long double h;
    long double ru [3], rv [3], ruu [3], ruv [3], rvv [3];
    long double nn [3];
    long double nnrm;
    long double un [3];
    long double E, F, G, L, M, N;
    long double den;
    long double kk;
    long double hh;
    int i;
    ml_surf_func_t fs [3];
    if (ML_UNLIKELY (!rx || !ry || !rz || !K || !H)) return ML_ERR_INVALID_ARG;
    if (ml_isnan (u) || ml_isnan (v) || !ml_isfinite (u) || !ml_isfinite (v)) {
        *K = ml_make_nan ();
        *H = ml_make_nan ();
        return ML_ERR_NAN_INPUT;
    }
    fs [0] = rx;
    fs [1] = ry;
    fs [2] = rz;
    /* Second-order stencils in double precision: h = 1e-4 balances
     * truncation h^2 against rounding eps/h^2. */
    h = 1e-4L * (1.0L + fabsl ((long double) u) + fabsl ((long double) v));
    for (i = 0; i < 3; i++) {
        double f0 = fs [i](u, v);
        double fp_u = fs [i](u + (double) h, v);
        double fm_u = fs [i](u - (double) h, v);
        double fp_v = fs [i](u, v + (double) h);
        double fm_v = fs [i](u, v - (double) h);
        double fpp = fs [i](u + (double) h, v + (double) h);
        double fpm = fs [i](u + (double) h, v - (double) h);
        double fmp = fs [i](u - (double) h, v + (double) h);
        double fmm = fs [i](u - (double) h, v - (double) h);
        if (ml_isnan (f0) || ml_isnan (fp_u) || ml_isnan (fm_u) || ml_isnan (fp_v) || ml_isnan (fm_v) ||
            ml_isnan (fpp) || ml_isnan (fpm) || ml_isnan (fmp) || ml_isnan (fmm)) {
            *K = ml_make_nan ();
            *H = ml_make_nan ();
            return ML_ERR_NAN_INPUT;
        }
        if (!ml_isfinite (f0) || !ml_isfinite (fp_u) || !ml_isfinite (fm_u) || !ml_isfinite (fp_v) || !ml_isfinite (fm_v)) {
            *K = ml_make_nan ();
            *H = ml_make_nan ();
            return ML_ERR_NAN_INPUT;
        }
        ru [i] = ((long double) fp_u - (long double) fm_u) / (2.0L * h);
        rv [i] = ((long double) fp_v - (long double) fm_v) / (2.0L * h);
        ruu [i] = ((long double) fp_u - 2.0L * (long double) f0 + (long double) fm_u) / (h * h);
        rvv [i] = ((long double) fp_v - 2.0L * (long double) f0 + (long double) fm_v) / (h * h);
        ruv [i] = ((long double) fpp - (long double) fpm - (long double) fmp + (long double) fmm) / (4.0L * h * h);
    }
    nn [0] = ru [1] * rv [2] - ru [2] * rv [1];
    nn [1] = ru [2] * rv [0] - ru [0] * rv [2];
    nn [2] = ru [0] * rv [1] - ru [1] * rv [0];
    nnrm = __builtin_sqrtl (nn [0] * nn [0] + nn [1] * nn [1] + nn [2] * nn [2]);
    if (!(nnrm > 0.0L) || !ml_isfinite ((double) nnrm)) {
        *K = ml_make_nan ();
        *H = ml_make_nan ();
        return ML_ERR_SINGULAR;
    }
    for (i = 0; i < 3; i++) un [i] = nn [i] / nnrm;
    E = ru [0] * ru [0] + ru [1] * ru [1] + ru [2] * ru [2];
    F = ru [0] * rv [0] + ru [1] * rv [1] + ru [2] * rv [2];
    G = rv [0] * rv [0] + rv [1] * rv [1] + rv [2] * rv [2];
    L = ruu [0] * un [0] + ruu [1] * un [1] + ruu [2] * un [2];
    M = ruv [0] * un [0] + ruv [1] * un [1] + ruv [2] * un [2];
    N = rvv [0] * un [0] + rvv [1] * un [1] + rvv [2] * un [2];
    den = E * G - F * F;
    if (!(den > 0.0L) || !ml_isfinite ((double) den)) {
        *K = ml_make_nan ();
        *H = ml_make_nan ();
        return ML_ERR_SINGULAR;
    }
    kk = (L * N - M * M) / den;
    hh = (L * G - 2.0L * M * F + N * E) / (2.0L * den);
    *K = ml_isfinite ((double) kk) ? (double) kk : ml_make_nan ();
    *H = ml_isfinite ((double) hh) ? (double) hh : ml_make_nan ();
    return ML_SUCCESS;
}
