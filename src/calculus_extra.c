#include "ml_compiler.h"
#include "ml_calculus_extra.h"
#include <stddef.h>

/* v12R2 calculus depth pack. Stateless, thread-safe, C99 only. */

#define CALC_MAXN 256
static double central_d (ml_func_t f, double x, double h) {
    return (f (x + h) - f (x - h)) / (2.0 * h);
}
ML_API double ml_derivative_richardson (ml_func_t f, double x, double h,
                                        double *err_out) {
    double d0, d1, d2, e0, e1, r;
    if (h <= 0.0) h = 1e-5;
    d0 = central_d (f, x, h);
    d1 = central_d (f, x, h * 0.5);
    d2 = central_d (f, x, h * 0.25);
    e0 = (4.0 * d1 - d0) / 3.0;
    e1 = (4.0 * d2 - d1) / 3.0;
    r = (16.0 * e1 - e0) / 15.0;
    if (err_out) *err_out = ml_fabs (r - e1);
    return r;
}
ML_API ml_status_t ml_gradient_rn (ml_rn_func_t f, const double *x, int n,
                                   double *grad, double h, void *ctx) {
    double xp [CALC_MAXN];
    int i;
    if (ML_UNLIKELY (!f || !x || !grad || n <= 0)) return ML_ERR_INVALID_ARG;
    if (ML_UNLIKELY (n > CALC_MAXN)) return ML_ERR_WORKSPACE;
    if (h <= 0.0) h = 1e-6;
    for (i = 0; i < n; i++) xp [i] = x [i];
    for (i = 0; i < n; i++) {
        double fp, fm, old = xp [i];
        xp [i] = old + h;
        fp = f (xp, n, ctx);
        xp [i] = old - h;
        fm = f (xp, n, ctx);
        xp [i] = old;
        grad [i] = (fp - fm) / (2.0 * h);
    }
    return ML_SUCCESS;
}
ML_API ml_status_t ml_jacobian_rn (ml_rn_vec_func_t f, const double *x, int n,
                                   int m, double *J, double h, void *ctx) {
    double xp [CALC_MAXN], fp [CALC_MAXN], fm [CALC_MAXN];
    int i, j;
    if (ML_UNLIKELY (!f || !x || !J || n <= 0 || m <= 0)) return ML_ERR_INVALID_ARG;
    if (ML_UNLIKELY (n > CALC_MAXN || m > CALC_MAXN)) return ML_ERR_WORKSPACE;
    if (h <= 0.0) h = 1e-6;
    for (i = 0; i < n; i++) xp [i] = x [i];
    for (j = 0; j < n; j++) {
        double old = xp [j];
        xp [j] = old + h;
        f (xp, n, fp, m, ctx);
        xp [j] = old - h;
        f (xp, n, fm, m, ctx);
        xp [j] = old;
        for (i = 0; i < m; i++)
            J [(size_t) i * (size_t) n + (size_t) j] = (fp [i] - fm [i]) / (2.0 * h);
    }
    return ML_SUCCESS;
}
ML_API ml_status_t ml_hessian_rn (ml_rn_func_t f, const double *x, int n,
                                  double *H, double h, void *ctx) {
    double xp [CALC_MAXN];
    int i, j;
    if (ML_UNLIKELY (!f || !x || !H || n <= 0)) return ML_ERR_INVALID_ARG;
    if (ML_UNLIKELY (n > CALC_MAXN)) return ML_ERR_WORKSPACE;
    if (h <= 0.0) h = 1e-4;
    for (i = 0; i < n; i++) xp [i] = x [i];
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            double xi = xp [i], xj = xp [j];
            if (i == j) {
                double fpp, f0, fmm;
                xp [i] = xi + h;
                fpp = f (xp, n, ctx);
                xp [i] = xi;
                f0 = f (xp, n, ctx);
                xp [i] = xi - h;
                fmm = f (xp, n, ctx);
                xp [i] = xi;
                H [(size_t) i * (size_t) n + (size_t) j] =
                    (fpp - 2.0 * f0 + fmm) / (h * h);
            } else {
                double fpp, fpm, fmp, fmm;
                xp [i] = xi + h;
                xp [j] = xj + h;
                fpp = f (xp, n, ctx);
                xp [i] = xi + h;
                xp [j] = xj - h;
                fpm = f (xp, n, ctx);
                xp [i] = xi - h;
                xp [j] = xj + h;
                fmp = f (xp, n, ctx);
                xp [i] = xi - h;
                xp [j] = xj - h;
                fmm = f (xp, n, ctx);
                xp [i] = xi;
                xp [j] = xj;
                H [(size_t) i * (size_t) n + (size_t) j] =
                    (fpp - fpm - fmp + fmm) / (4.0 * h * h);
            }
        }
    }
    return ML_SUCCESS;
}
ML_API double ml_directional_derivative (ml_rn_func_t f, const double *x, int n,
                                         const double *dir, double h, void *ctx) {
    double xp [CALC_MAXN], xm [CALC_MAXN];
    int i;
    if (ML_UNLIKELY (!f || !x || !dir || n <= 0 || n > CALC_MAXN)) return ml_make_nan ();
    if (h <= 0.0) h = 1e-6;
    for (i = 0; i < n; i++) {
        xp [i] = x [i] + h * dir [i];
        xm [i] = x [i] - h * dir [i];
    }
    return (f (xp, n, ctx) - f (xm, n, ctx)) / (2.0 * h);
}
ML_API double ml_limit (ml_func_t f, double x0, int side, double tol,
                        double *err_out) {
    double s [5];
    int k;
    if (ML_UNLIKELY (!f)) return ml_make_nan ();
    (void) tol;
    if (side == 0) {
        double l = ml_limit (f, x0, -1, tol, err_out);
        double r = ml_limit (f, x0, +1, tol, err_out);
        if (err_out) *err_out = ml_fabs (l - r) * 0.5;
        return 0.5 * (l + r);
    }
    for (k = 0; k < 5; k++) {
        double h = 1e-2;
        int j;
        for (j = 0; j < k; j++) h *= 0.5;
        s [k] = f (x0 + (double) side * h);
    }
    /* Richardson for geometric p = 1/2, applied twice, then Aitken. */
    {
        double t [4], u [3];
        for (k = 1; k < 5; k++) t [k - 1] = 2.0 * s [k] - s [k - 1];
        for (k = 1; k < 4; k++) u [k - 1] = 2.0 * t [k] - t [k - 1];
        if (err_out) *err_out = ml_fabs (u [2] - u [1]);
        return ml_aitken (u, 3);
    }
}
ML_API double ml_aitken (const double *s, int n) {
    double e [64][64];
    int i, j, k;
    if (ML_UNLIKELY (!s || n <= 0)) return ml_make_nan ();
    if (n > 64) n = 64;
    if (n == 1) return s [0];
    for (i = 0; i < n; i++) e [i][0] = s [i];
    for (j = 1; j < n; j++) {
        for (k = 0; k + j < n; k++) {
            double a = (j == 1) ? 0.0 : e [k + 1][j - 2];
            double denom = e [k + 1][j - 1] - e [k][j - 1];
            e [k][j] = (denom == 0.0) ? 1e300 : a + 1.0 / denom;
        }
    }
    j = 2 * ((n - 1) / 2);
    return e [0][j];
}
ML_API double ml_richardson_geom (const double *s, int n, double q) {
    if (ML_UNLIKELY (!s || n < 2)) return (s && n > 0) ? s [n - 1] : ml_make_nan ();
    if (q == 1.0) return s [n - 1];
    return (s [n - 1] - q * s [n - 2]) / (1.0 - q);
}
ML_API ml_status_t ml_pade (const double *c, int L, int M, double *p, double *q) {
    double A [40][41];
    int i, j, k;
    if (ML_UNLIKELY (!c || L < 0 || M < 0 || (L > 0 && !p) || !q)) return ML_ERR_INVALID_ARG;
    if (ML_UNLIKELY (L + M > 39)) return ML_ERR_WORKSPACE;
    if (L > 0 && p) p [0] = 0.0;
    q [0] = 1.0;
    if (M == 0) {
        for (k = 0; k <= L; k++) p [k] = c [k];
        return ML_SUCCESS;
    }
    /* Build M x (M+1) augmented system for q[1..M]. */
    for (i = 0; i < M; i++) {
        int row = L + 1 + i;
        for (j = 1; j <= M; j++) {
            int idx = row - j;
            A [i][j - 1] = (idx >= 0) ? c [idx] : 0.0;
        }
        A [i][M] = (row >= 0) ? -c [row] : 0.0;
    }
    for (j = 0; j < M; j++) {
        int piv = j;
        double best = ml_fabs (A [j][j]);
        for (i = j + 1; i < M; i++) {
            double v = ml_fabs (A [i][j]);
            if (v > best) {
                best = v;
                piv = i;
            }
        }
        if (best == 0.0) return ML_ERR_SINGULAR;
        if (piv != j)
            for (k = j; k <= M; k++) {
                double t = A [j][k];
                A [j][k] = A [piv][k];
                A [piv][k] = t;
            }
        {
            double d = A [j][j];
            for (k = j; k <= M; k++) A [j][k] /= d;
        }
        for (i = 0; i < M; i++) {
            double fac;
            if (i == j) continue;
            fac = A [i][j];
            if (fac == 0.0) continue;
            for (k = j; k <= M; k++) A [i][k] -= fac * A [j][k];
        }
    }
    for (j = 1; j <= M; j++) q [j] = A [j - 1][M];
    for (k = 0; k <= L; k++) {
        double s = 0.0;
        int jmax = (k < M) ? k : M;
        for (j = 0; j <= jmax; j++) {
            int idx = k - j;
            if (idx >= 0) s += c [idx] * q [j];
        }
        p [k] = s;
    }
    return ML_SUCCESS;
}
ML_API double ml_pade_eval (const double *p, const double *q, int L, int M,
                            double x) {
    double num = 0.0, den = 0.0;
    int k;
    if (ML_UNLIKELY (!p || !q)) return ml_make_nan ();
    for (k = L; k >= 0; k--) num = num * x + p [k];
    for (k = M; k >= 0; k--) den = den * x + q [k];
    if (den == 0.0) return ml_make_nan ();
    return num / den;
}
/* -------------------------------------------------------------- integrals - */

static double simp_w (int i, int n) {
    if (i == 0 || i == n) return 1.0;
    return (i & 1) ? 4.0 : 2.0;
}
ML_API double ml_integral_2d (ml_r2_func_t f, double xa, double xb, double ya,
                              double yb, int nx, int ny, void *ctx) {
    double hx, hy, s = 0.0;
    int i, j;
    if (ML_UNLIKELY (!f)) return ml_make_nan ();
    if (nx < 2) nx = 2;
    if (ny < 2) ny = 2;
    if (nx & 1) nx++;
    if (ny & 1) ny++;
    hx = (xb - xa) / (double) nx;
    hy = (yb - ya) / (double) ny;
    for (i = 0; i <= nx; i++) {
        double xi = xa + (double) i * hx, wx = simp_w (i, nx);
        for (j = 0; j <= ny; j++) {
            double yj = ya + (double) j * hy, wy = simp_w (j, ny);
            s += wx * wy * f (xi, yj, ctx);
        }
    }
    return s * hx * hy / 9.0;
}
ML_API double ml_integral_2d_var (ml_r2_func_t f, double xa, double xb,
                                  ml_r1c_func_t ylo, ml_r1c_func_t yhi,
                                  int nx, int ny, void *ctx) {
    double hx, s = 0.0;
    int i, j;
    if (ML_UNLIKELY (!f || !ylo || !yhi)) return ml_make_nan ();
    if (nx < 2) nx = 2;
    if (ny < 2) ny = 2;
    if (nx & 1) nx++;
    if (ny & 1) ny++;
    hx = (xb - xa) / (double) nx;
    for (i = 0; i <= nx; i++) {
        double xi = xa + (double) i * hx, wx = simp_w (i, nx);
        double a = ylo (xi, ctx), b = yhi (xi, ctx), hy = (b - a) / (double) ny, inner = 0.0;
        for (j = 0; j <= ny; j++)
            inner += simp_w (j, ny) * f (xi, a + (double) j * hy, ctx);
        s += wx * (hy / 3.0) * inner;
    }
    return (hx / 3.0) * s;
}
ML_API double ml_integral_3d (ml_r3c_func_t f, double xa, double xb, double ya,
                              double yb, double za, double zb, int n, void *ctx) {
    double hx, hy, hz, s = 0.0;
    int i, j, k;
    if (ML_UNLIKELY (!f)) return ml_make_nan ();
    if (n < 2) n = 2;
    if (n & 1) n++;
    hx = (xb - xa) / (double) n;
    hy = (yb - ya) / (double) n;
    hz = (zb - za) / (double) n;
    for (i = 0; i <= n; i++) {
        double xi = xa + (double) i * hx, wx = simp_w (i, n);
        for (j = 0; j <= n; j++) {
            double yj = ya + (double) j * hy, wy = simp_w (j, n);
            for (k = 0; k <= n; k++) {
                double zk = za + (double) k * hz, wz = simp_w (k, n);
                s += wx * wy * wz * f (xi, yj, zk, ctx);
            }
        }
    }
    return s * hx * hy * hz / 27.0;
}
/* ------------------------------------------------------------------ curves - */

static void curve_pt (ml_func_t rx, ml_func_t ry, ml_func_t rz, double t,
                      double r [3]) {
    r [0] = rx (t);
    r [1] = ry (t);
    r [2] = rz (t);
}
static void curve_d (ml_func_t rx, ml_func_t ry, ml_func_t rz, double t, double h,
                     double d [3]) {
    d [0] = (rx (t + h) - rx (t - h)) / (2.0 * h);
    d [1] = (ry (t + h) - ry (t - h)) / (2.0 * h);
    d [2] = (rz (t + h) - rz (t - h)) / (2.0 * h);
}
ML_API double ml_arc_length (ml_func_t rx, ml_func_t ry, ml_func_t rz,
                             double a, double b, int n) {
    double h, s = 0.0;
    int i;
    if (ML_UNLIKELY (!rx || !ry || !rz)) return ml_make_nan ();
    if (n < 2) n = 2;
    if (n & 1) n++;
    h = (b - a) / (double) n;
    for (i = 0; i <= n; i++) {
        double d [3], t = a + (double) i * h;
        curve_d (rx, ry, rz, t, 1e-6, d);
        s += simp_w (i, n) * ml_sqrt (d [0] * d [0] + d [1] * d [1] + d [2] * d [2]);
    }
    return s * h / 3.0;
}
ML_API double ml_line_integral_scalar (ml_r3c_func_t f, ml_func_t rx,
                                       ml_func_t ry, ml_func_t rz, double a,
                                       double b, int n, void *ctx) {
    double h, s = 0.0;
    int i;
    if (ML_UNLIKELY (!f || !rx || !ry || !rz)) return ml_make_nan ();
    if (n < 2) n = 2;
    if (n & 1) n++;
    h = (b - a) / (double) n;
    for (i = 0; i <= n; i++) {
        double r [3], d [3], t = a + (double) i * h, dl;
        curve_pt (rx, ry, rz, t, r);
        curve_d (rx, ry, rz, t, 1e-6, d);
        dl = ml_sqrt (d [0] * d [0] + d [1] * d [1] + d [2] * d [2]);
        s += simp_w (i, n) * f (r [0], r [1], r [2], ctx) * dl;
    }
    return s * h / 3.0;
}
ML_API double ml_line_integral_vector (ml_r3c_func_t fx, ml_r3c_func_t fy,
                                       ml_r3c_func_t fz, ml_func_t rx,
                                       ml_func_t ry, ml_func_t rz, double a,
                                       double b, int n, void *ctx) {
    double h, s = 0.0;
    int i;
    if (ML_UNLIKELY (!fx || !fy || !fz || !rx || !ry || !rz)) return ml_make_nan ();
    if (n < 2) n = 2;
    if (n & 1) n++;
    h = (b - a) / (double) n;
    for (i = 0; i <= n; i++) {
        double r [3], d [3], t = a + (double) i * h;
        curve_pt (rx, ry, rz, t, r);
        curve_d (rx, ry, rz, t, 1e-6, d);
        s += simp_w (i, n) * (fx (r [0], r [1], r [2], ctx) * d [0] + fy (r [0], r [1], r [2], ctx) * d [1] + fz (r [0], r [1], r [2], ctx) * d [2]);
    }
    return s * h / 3.0;
}
/* ----------------------------------------------------------------- surfaces */

static double surf_elem (ml_uv_func_t rx, ml_uv_func_t ry, ml_uv_func_t rz, double u,
                         double v) {
    double h = 1e-6;
    double ru [3], rv [3], cx, cy, cz;
    ru [0] = (rx (u + h, v) - rx (u - h, v)) / (2.0 * h);
    ru [1] = (ry (u + h, v) - ry (u - h, v)) / (2.0 * h);
    ru [2] = (rz (u + h, v) - rz (u - h, v)) / (2.0 * h);
    rv [0] = (rx (u, v + h) - rx (u, v - h)) / (2.0 * h);
    rv [1] = (ry (u, v + h) - ry (u, v - h)) / (2.0 * h);
    rv [2] = (rz (u, v + h) - rz (u, v - h)) / (2.0 * h);
    cx = ru [1] * rv [2] - ru [2] * rv [1];
    cy = ru [2] * rv [0] - ru [0] * rv [2];
    cz = ru [0] * rv [1] - ru [1] * rv [0];
    return ml_sqrt (cx * cx + cy * cy + cz * cz);
}
ML_API double ml_surface_integral_scalar (ml_r3c_func_t f, ml_uv_func_t rx,
                                          ml_uv_func_t ry, ml_uv_func_t rz, double ua,
                                          double ub, double va, double vb,
                                          int nu, int nv, void *ctx) {
    double hu, hv, s = 0.0;
    int i, j;
    if (ML_UNLIKELY (!f || !rx || !ry || !rz)) return ml_make_nan ();
    if (nu < 2) nu = 2;
    if (nv < 2) nv = 2;
    if (nu & 1) nu++;
    if (nv & 1) nv++;
    hu = (ub - ua) / (double) nu;
    hv = (vb - va) / (double) nv;
    for (i = 0; i <= nu; i++) {
        double u = ua + (double) i * hu, wu = simp_w (i, nu);
        for (j = 0; j <= nv; j++) {
            double v = va + (double) j * hv, wv = simp_w (j, nv);
            s += wu * wv * f (rx (u, v), ry (u, v), rz (u, v), ctx) * surf_elem (rx, ry, rz, u, v);
        }
    }
    return s * hu * hv / 9.0;
}
