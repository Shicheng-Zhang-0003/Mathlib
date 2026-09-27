#include "ml_compiler.h"
#include "ml_calculus.h"
#include "ml_orthogonal.h"

ML_API ml_status_t ml_gauss_legendre(int n, double *xs, double *ws) {
    if (ML_UNLIKELY(!xs || !ws || n <= 0 || n > 1024)) return ML_ERR_INVALID_ARG;
    for (int i = 1; i <= n; i++) {
        /* Tricomi initial guess: cos(pi*(i-0.25)/(n+0.5)) */
        long double x = __builtin_cosl(3.14159265358979323846264338327950288L
                        * ((long double)i - 0.25L) / ((long double)n + 0.5L));
        for (int it = 0; it < 50; it++) {
            long double p0 = 1.0L, p1 = x, p2 = 0.0L;
            for (int k = 1; k < n; k++) {
                p2 = ((2.0L * (long double)k + 1.0L) * x * p1
                      - (long double)k * p0) / (long double)(k + 1);
                p0 = p1; p1 = p2;
            }
            /* deriv via (x^2-1)P' = n(xP_n - P_{n-1}) */
            long double dp = (long double)n * (x * p1 - p0) / (x * x - 1.0L);
            long double dx = p1 / dp;
            x -= dx;
            if (x > 1.0L) x = 1.0L;
            if (x < -1.0L) x = -1.0L;
            if (__builtin_fabsl(dx) < 1e-20L) break;
        }
        {
            long double p0 = 1.0L, p1 = x, p2 = 0.0L;
            for (int k = 1; k < n; k++) {
                p2 = ((2.0L * (long double)k + 1.0L) * x * p1
                      - (long double)k * p0) / (long double)(k + 1);
                p0 = p1; p1 = p2;
            }
            /* Residual check: Newton must have driven P_n(x) to ~0.
             * ml_types.h defines no ML_ERR_CONVERGENCE, so a residual
             * above 1e-14 is reported as ML_ERR_INTERNAL. */
            if (__builtin_fabsl(p1) > 1e-14L) {
                return ML_ERR_INTERNAL;
            }
            long double dp = (long double)n * (x * p1 - p0) / (x * x - 1.0L);
            long double w = 2.0L / ((1.0L - x * x) * dp * dp);
            xs[i - 1] = (double)x;
            ws[i - 1] = (double)w;
            if (!ml_isfinite(xs[i-1]) || !ml_isfinite(ws[i-1])) return ML_ERR_INTERNAL;
        }
    }
    /* Sort ascending (Newton fills descending for this guess order). */
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (xs[j] < xs[i]) {
                double t = xs[i]; xs[i] = xs[j]; xs[j] = t;
                t = ws[i]; ws[i] = ws[j]; ws[j] = t;
            }
        }
    }
    return ML_SUCCESS;
}

ML_API double ml_gauss_legendre_integral(ml_func_t f, double a, double b, int n) {
    if (ML_UNLIKELY(!f || n <= 0 || n > 1024)) return ml_make_nan();
    if (ML_UNLIKELY(ml_isnan(a) || ml_isnan(b) || ml_isinf(a) || ml_isinf(b))) return ml_make_nan();
    if (a == b) return 0.0;
    /* Signed integral: a>b swaps limits and negates (instead of NaN). */
    int negate = 0;
    if (a > b) {
        double t = a; a = b; b = t;
        negate = 1;
    }
    {
        double xs[1024], ws[1024];
        if (ml_gauss_legendre(n, xs, ws) != ML_SUCCESS) return ml_make_nan();
        long double c = ((long double)b + (long double)a) * 0.5L;
        long double h = ((long double)b - (long double)a) * 0.5L;
        long double s = 0.0L;
        for (int i = 0; i < n; i++) {
            double xi = (double)(c + h * (long double)xs[i]);
            double fi = f(xi);
            if (!ml_isfinite(fi)) return ml_make_nan();
            s += (long double)fi * (long double)ws[i];
        }
        double r = (double)(s * h);
        if (!ml_isfinite(r)) return ml_make_nan();
        return (negate != 0) ? -r : r;
    }
}

ML_API double ml_lagrange_interp(const double *x, const double *y, int n, double x0) {
    if (ML_UNLIKELY(!x || !y || n <= 0 || n > 64)) return ml_make_nan();
    if (ML_UNLIKELY(ml_isnan(x0) || ml_isinf(x0))) return ml_make_nan();
    for (int i = 0; i < n; i++) {
        if (ML_UNLIKELY(!ml_isfinite(x[i]) || !ml_isfinite(y[i]))) return ml_make_nan();
        if (x0 == x[i]) return y[i];
        for (int j = i + 1; j < n; j++) {
            if (x[i] == x[j]) return ml_make_nan();
        }
    }
    {
        /* Barycentric form: w_i = 1/prod_{j!=i}(x_i-x_j), long double. */
        long double w[64];
        for (int i = 0; i < n; i++) {
            long double d = 1.0L;
            for (int j = 0; j < n; j++) {
                if (j == i) continue;
                d *= (long double)x[i] - (long double)x[j];
            }
            if (d == 0.0L) return ml_make_nan();
            w[i] = 1.0L / d;
        }
        long double num = 0.0L, den = 0.0L;
        for (int i = 0; i < n; i++) {
            long double t = w[i] / ((long double)x0 - (long double)x[i]);
            num += t * (long double)y[i];
            den += t;
        }
        if (den == 0.0L) return ml_make_nan();
        double r = (double)(num / den);
        return ml_isfinite(r) ? r : ml_make_nan();
    }
}

ML_API ml_status_t ml_cubic_spline_natural(const double *x, const double *y, int n, double *m2) {
    if (ML_UNLIKELY(!x || !y || !m2 || n < 2 || n > 1024)) return ML_ERR_INVALID_ARG;
    for (int i = 0; i < n; i++) {
        if (ML_UNLIKELY(!ml_isfinite(x[i]) || !ml_isfinite(y[i]))) return ML_ERR_NAN_INPUT;
        if (i > 0 && !(x[i] > x[i-1])) return ML_ERR_INVALID_ARG;
    }
    m2[0] = 0.0; m2[n-1] = 0.0;
    if (n == 2) return ML_SUCCESS;
    {
        /* Thomas tridiagonal for natural spline second derivatives.
         * Stack-local scratch (thread-safe; no static storage). */
        double h[1024], rhs[1024], cp[1024], dp2[1024];
        for (int i = 0; i < n - 1; i++) h[i] = x[i+1] - x[i];
        for (int i = 1; i < n - 1; i++) {
            rhs[i] = 6.0 * ((y[i+1] - y[i]) / h[i] - (y[i] - y[i-1]) / h[i-1]);
        }
        /* Forward: 2(h_{i-1}+h_i) m_i + h_i m_{i+1} ... */
        cp[1] = h[1] / (2.0 * (h[0] + h[1]));
        dp2[1] = rhs[1] / (2.0 * (h[0] + h[1]));
        for (int i = 2; i < n - 1; i++) {
            double den = 2.0 * (h[i-1] + h[i]) - h[i-1] * cp[i-1];
            if (den == 0.0) return ML_ERR_SINGULAR;
            cp[i] = h[i] / den;
            dp2[i] = (rhs[i] - h[i-1] * dp2[i-1]) / den;
        }
        m2[n-2] = dp2[n-2];
        for (int i = n - 3; i >= 1; i--) m2[i] = dp2[i] - cp[i] * m2[i+1];
        for (int i = 1; i < n - 1; i++) {
            if (!ml_isfinite(m2[i])) return ML_ERR_INTERNAL;
        }
        return ML_SUCCESS;
    }
}

ML_API double ml_cubic_spline_eval(const double *x, const double *y, const double *m2, int n, double x0) {
    if (ML_UNLIKELY(!x || !y || !m2 || n < 2)) return ml_make_nan();
    if (ML_UNLIKELY(ml_isnan(x0) || ml_isinf(x0))) return ml_make_nan();
    if (x0 < x[0] || x0 > x[n-1]) return ml_make_nan();
    {
        int k = 0;
        while (k + 1 < n - 1 && x0 > x[k+1]) k++;
        long double h = (long double)x[k+1] - (long double)x[k];
        if (h <= 0.0L) return ml_make_nan();
        long double a = ((long double)x[k+1] - (long double)x0) / h;
        long double b = ((long double)x0 - (long double)x[k]) / h;
        long double s = a * (long double)y[k] + b * (long double)y[k+1]
            + ((a*a*a - a) * (long double)m2[k] + (b*b*b - b) * (long double)m2[k+1])
              * h * h / 6.0L;
        double r = (double)s;
        return ml_isfinite(r) ? r : ml_make_nan();
    }
}

ML_API void ml_grad3(ml_vec3_func_t f3[3], double x, double y, double z, double h,
                     double *gx, double *gy, double *gz) {
    if (!gx || !gy || !gz) return;
    *gx = *gy = *gz = ml_make_nan();
    if (!f3 || !f3[0] || !f3[1] || !f3[2]) return;
    if (!ml_isfinite(x) || !ml_isfinite(y) || !ml_isfinite(z)) return;
    /* Standardized per-axis steps: caller h>0 is used for all axes,
     * otherwise sqrt(eps)*(1+|coord|) per axis. */
    double hx, hy, hz;
    if ((h > 0.0) && ml_isfinite(h)) {
        hx = h; hy = h; hz = h;
    } else {
        const double eps = 1.4901161193847656e-08;
        hx = eps * (1.0 + ml_fabs(x));
        hy = eps * (1.0 + ml_fabs(y));
        hz = eps * (1.0 + ml_fabs(z));
    }
    /* True partials of F:R^3->R^3 along each component's own axis:
     * gx=dF0/dx, gy=dF1/dy, gz=dF2/dz. */
    {
        double f0p = f3[0](x + hx, y, z);
        double f0m = f3[0](x - hx, y, z);
        double f1p = f3[1](x, y + hy, z);
        double f1m = f3[1](x, y - hy, z);
        double f2p = f3[2](x, y, z + hz);
        double f2m = f3[2](x, y, z - hz);
        if (!ml_isfinite(f0p) || !ml_isfinite(f0m) ||
            !ml_isfinite(f1p) || !ml_isfinite(f1m) ||
            !ml_isfinite(f2p) || !ml_isfinite(f2m)) {
            return;
        }
        *gx = (f0p - f0m) / (2.0 * hx);
        *gy = (f1p - f1m) / (2.0 * hy);
        *gz = (f2p - f2m) / (2.0 * hz);
    }
}

ML_API double ml_div3(ml_vec3_func_t f3[3], double x, double y, double z, double h) {
    double gx, gy, gz;
    ml_grad3(f3, x, y, z, h, &gx, &gy, &gz);
    if (!ml_isfinite(gx) || !ml_isfinite(gy) || !ml_isfinite(gz)) return ml_make_nan();
    return gx + gy + gz;
}

ML_API void ml_curl3(ml_vec3_func_t f3[3], double x, double y, double z, double h,
                     double *cx, double *cy, double *cz) {
    if (!cx || !cy || !cz) return;
    *cx = *cy = *cz = ml_make_nan();
    if (!f3 || !f3[0] || !f3[1] || !f3[2]) return;
    if (!ml_isfinite(x) || !ml_isfinite(y) || !ml_isfinite(z)) return;
    /* Standardized per-axis steps (same rule as ml_grad3). */
    double hx, hy, hz;
    if ((h > 0.0) && ml_isfinite(h)) {
        hx = h; hy = h; hz = h;
    } else {
        const double eps = 1.4901161193847656e-08;
        hx = eps * (1.0 + ml_fabs(x));
        hy = eps * (1.0 + ml_fabs(y));
        hz = eps * (1.0 + ml_fabs(z));
    }
    /* Curl via true off-axis partials of F(x,y,z):
     * cx=dF2/dy-dF1/dz, cy=dF0/dz-dF2/dx, cz=dF1/dx-dF0/dy. */
    {
        double f2_yp = f3[2](x, y + hy, z);
        double f2_ym = f3[2](x, y - hy, z);
        double f1_zp = f3[1](x, y, z + hz);
        double f1_zm = f3[1](x, y, z - hz);
        double f0_zp = f3[0](x, y, z + hz);
        double f0_zm = f3[0](x, y, z - hz);
        double f2_xp = f3[2](x + hx, y, z);
        double f2_xm = f3[2](x - hx, y, z);
        double f1_xp = f3[1](x + hx, y, z);
        double f1_xm = f3[1](x - hx, y, z);
        double f0_yp = f3[0](x, y + hy, z);
        double f0_ym = f3[0](x, y - hy, z);
        if (!ml_isfinite(f2_yp) || !ml_isfinite(f2_ym) ||
            !ml_isfinite(f1_zp) || !ml_isfinite(f1_zm) ||
            !ml_isfinite(f0_zp) || !ml_isfinite(f0_zm) ||
            !ml_isfinite(f2_xp) || !ml_isfinite(f2_xm) ||
            !ml_isfinite(f1_xp) || !ml_isfinite(f1_xm) ||
            !ml_isfinite(f0_yp) || !ml_isfinite(f0_ym)) {
            return;
        }
        {
            double dwy = (f2_yp - f2_ym) / (2.0 * hy);
            double dvz = (f1_zp - f1_zm) / (2.0 * hz);
            double duz = (f0_zp - f0_zm) / (2.0 * hz);
            double dwx = (f2_xp - f2_xm) / (2.0 * hx);
            double dvx = (f1_xp - f1_xm) / (2.0 * hx);
            double duy = (f0_yp - f0_ym) / (2.0 * hy);
            if (!ml_isfinite(dwy) || !ml_isfinite(dvz) || !ml_isfinite(duz) ||
                !ml_isfinite(dwx) || !ml_isfinite(dvx) || !ml_isfinite(duy)) {
                return;
            }
            *cx = dwy - dvz;
            *cy = duz - dwx;
            *cz = dvx - duy;
        }
    }
}
