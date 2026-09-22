#include "ml_compiler.h"
#include "ml_calculus.h"
#include "ml_orthogonal.h"

ML_API ml_status_t ml_gauss_legendre(int n, double *xs, double *ws) {
    if (ML_UNLIKELY(!xs || !ws || n <= 0 || n > 64)) return ML_ERR_INVALID_ARG;
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
    if (ML_UNLIKELY(!f || n <= 0 || n > 64)) return ml_make_nan();
    if (ML_UNLIKELY(ml_isnan(a) || ml_isnan(b) || ml_isinf(a) || ml_isinf(b))) return ml_make_nan();
    if (a == b) return 0.0;
    if (a > b) return ml_make_nan();
    {
        double xs[64], ws[64];
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
        return ml_isfinite(r) ? r : ml_make_nan();
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
        /* Thomas tridiagonal for natural spline second derivatives. */
        static double h[1024], rhs[1024], cp[1024], dp2[1024];
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

static double ml_fd1(ml_func_t f, double v, double h) {
    double fp = f(v + h), fm = f(v - h);
    return (fp - fm) / (2.0 * h);
}

ML_API void ml_grad3(ml_func_t f3[3], double x, double y, double z, double h,
                     double *gx, double *gy, double *gz) {
    (void)z;
    if (!gx || !gy || !gz) return;
    *gx = *gy = *gz = ml_make_nan();
    if (!f3 || !f3[0] || !f3[1] || !f3[2]) return;
    if (!ml_isfinite(x) || !ml_isfinite(y) || !ml_isfinite(z)) return;
    if (!(h > 0.0) || !ml_isfinite(h)) {
        h = 1.4901161193847656e-08 * (1.0 + ml_fabs(x) + ml_fabs(y) + ml_fabs(z));
    }
    /* Scalar-potential gradient would use one f; here F:R^3->R^3
     * Jacobian diagonal (du/dx, dv/dy, dw/dz) as grad-like check. */
    {
        double ax = ml_fabs(x), hxv = 1.4901161193847656e-08 * (1.0 + ax);
        double ay = ml_fabs(y), hyv = 1.4901161193847656e-08 * (1.0 + ay);
        double az = ml_fabs(z), hzv = 1.4901161193847656e-08 * (1.0 + az);
        (void)h;
        /* Central differences of each component along its own axis. */
        double u0 = f3[0](x), u1 = 0.0, v0 = 0.0, v1 = 0.0, w0 = 0.0, w1 = 0.0;
        (void)u0; (void)u1; (void)v0; (void)v1; (void)w0; (void)w1;
        /* Generic path: caller passes closures capturing y,z etc. via
         * file-static state; we evaluate along x/y/z shifts of each. */
        *gx = ml_fd1(f3[0], x, hxv);
        *gy = ml_fd1(f3[1], y, hyv);
        *gz = ml_fd1(f3[2], z, hzv);
    }
}

ML_API double ml_div3(ml_func_t f3[3], double x, double y, double z, double h) {
    double gx, gy, gz;
    ml_grad3(f3, x, y, z, h, &gx, &gy, &gz);
    if (!ml_isfinite(gx) || !ml_isfinite(gy) || !ml_isfinite(gz)) return ml_make_nan();
    return gx + gy + gz;
}

ML_API void ml_curl3(ml_func_t f3[3], double x, double y, double z, double h,
                     double *cx, double *cy, double *cz) {
    if (!cx || !cy || !cz) return;
    *cx = *cy = *cz = ml_make_nan();
    if (!f3 || !f3[0] || !f3[1] || !f3[2]) return;
    if (!ml_isfinite(x) || !ml_isfinite(y) || !ml_isfinite(z)) return;
    (void)h;
    {
        double ax = ml_fabs(x), ay = ml_fabs(y), az = ml_fabs(z);
        double hx = 1.4901161193847656e-08 * (1.0 + ax + ay + az);
        double hy = hx, hz = hx;
        /* Curl via central differences of off-axis partials. The func
         * array entries are 1-D slices provided by the caller. */
        double dwy = ml_fd1(f3[2], y, hy);
        double dvz = ml_fd1(f3[1], z, hz);
        double duz = ml_fd1(f3[0], z, hz);
        double dwx = ml_fd1(f3[2], x, hx);
        double dvx = ml_fd1(f3[1], x, hx);
        double duy = ml_fd1(f3[0], y, hy);
        if (!ml_isfinite(dwy) || !ml_isfinite(dvz) || !ml_isfinite(duz) ||
            !ml_isfinite(dwx) || !ml_isfinite(dvx) || !ml_isfinite(duy)) return;
        *cx = dwy - dvz;
        *cy = duz - dwx;
        *cz = dvx - duy;
    }
}
