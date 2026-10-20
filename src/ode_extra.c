#include "ml_compiler.h"
#include "ml_ode_extra.h"
/* ============================================================================
 * v12R2 ODE / BVP depth pack.
 * ========================================================================== */

static ml_status_t oe_solve (double *A, double *b, double *x, int n) {
    int i, j, k;
    for (k = 0; k < n; k++) {
        int piv = k;
        double mx = ml_fabs (A [(size_t) k * (size_t) n + (size_t) k]);
        for (i = k + 1; i < n; i++) {
            double v = ml_fabs (A [(size_t) i * (size_t) n + (size_t) k]);
            if (v > mx) {
                mx = v;
                piv = i;
            }
        }
        if (mx <= 0.0) return ML_ERR_SINGULAR;
        if (piv != k) {
            for (j = 0; j < n; j++) {
                double t = A [(size_t) k * (size_t) n + (size_t) j];
                A [(size_t) k * (size_t) n + (size_t) j] = A [(size_t) piv * (size_t) n + (size_t) j];
                A [(size_t) piv * (size_t) n + (size_t) j] = t;
            }
            {
                double t = b [k];
                b [k] = b [piv];
                b [piv] = t;
            }
        }
        {
            double d = A [(size_t) k * (size_t) n + (size_t) k];
            for (i = k + 1; i < n; i++) {
                double fct = A [(size_t) i * (size_t) n + (size_t) k] / d;
                b [i] -= fct * b [k];
                for (j = k; j < n; j++)
                    A [(size_t) i * (size_t) n + (size_t) j] -=
                        fct * A [(size_t) k * (size_t) n + (size_t) j];
            }
        }
    }
    for (i = n - 1; i >= 0; i--) {
        double s = b [i];
        for (j = i + 1; j < n; j++)
            s -= A [(size_t) i * (size_t) n + (size_t) j] * x [j];
        x [i] = s / A [(size_t) i * (size_t) n + (size_t) i];
    }
    return ML_SUCCESS;
}
ML_API ml_status_t ml_ode_rk4_sys (ml_sys_func_t f, double t0,
                                   const double *y0, int n, double dt, int steps,
                                   void *ctx, double *y_out) {
    double y [ML_ODE_MAXN], k1 [ML_ODE_MAXN], k2 [ML_ODE_MAXN], k3 [ML_ODE_MAXN],
        k4 [ML_ODE_MAXN], tmp [ML_ODE_MAXN];
    int i, s;
    if (ML_UNLIKELY (!f || !y0 || !y_out || n <= 0 || n > ML_ODE_MAXN))
        return ML_ERR_INVALID_ARG;
    for (i = 0; i < n; i++) y [i] = y0 [i];
    for (s = 0; s < steps; s++) {
        double t = t0 + dt * (double) s, h2 = 0.5 * dt;
        f (t, y, k1, n, ctx);
        for (i = 0; i < n; i++) tmp [i] = y [i] + h2 * k1 [i];
        f (t + h2, tmp, k2, n, ctx);
        for (i = 0; i < n; i++) tmp [i] = y [i] + h2 * k2 [i];
        f (t + h2, tmp, k3, n, ctx);
        for (i = 0; i < n; i++) tmp [i] = y [i] + dt * k3 [i];
        f (t + dt, tmp, k4, n, ctx);
        for (i = 0; i < n; i++)
            y [i] += dt / 6.0 * (k1 [i] + 2.0 * k2 [i] + 2.0 * k3 [i] + k4 [i]);
    }
    for (i = 0; i < n; i++) y_out [i] = y [i];
    return ML_SUCCESS;
}
/* Solve the implicit equation G(y)=0 with a finite-difference Newton, where
 * G(y) = y - base - theta*dt*f(t+dt,y). */
static ml_status_t oe_newton_step (ml_sys_func_t f, double tn, const double *base,
                                   int n, double theta_dt, void *ctx,
                                   double *y) {
    double g [ML_ODE_MAXN], fp [ML_ODE_MAXN], J [ML_ODE_MAXN * ML_ODE_MAXN],
        dy [ML_ODE_MAXN], rhs [ML_ODE_MAXN], yp [ML_ODE_MAXN];
    int i, j, it;
    for (it = 0; it < 50; it++) {
        double gn = 0.0;
        f (tn, y, fp, n, ctx);
        for (i = 0; i < n; i++) g [i] = y [i] - base [i] - theta_dt * fp [i];
        for (i = 0; i < n; i++) gn += g [i] * g [i];
        if (ml_sqrt (gn) < 1e-12) return ML_SUCCESS;
        for (j = 0; j < n; j++) {
            double yj = y [j], h = 1e-7 * (ml_fabs (yj) + 1.0);
            double fpp [ML_ODE_MAXN];
            for (i = 0; i < n; i++) yp [i] = y [i];
            yp [j] = yj + h;
            f (tn, yp, fpp, n, ctx);
            for (i = 0; i < n; i++)
                J [(size_t) i * (size_t) n + (size_t) j] =
                    ((i == j) ? 1.0 : 0.0) - theta_dt * (fpp [i] - fp [i]) / h;
        }
        for (i = 0; i < n; i++) rhs [i] = -g [i];
        if (oe_solve (J, rhs, dy, n) != ML_SUCCESS) return ML_ERR_SINGULAR;
        for (i = 0; i < n; i++) y [i] += dy [i];
    }
    return ML_SUCCESS;
}
ML_API ml_status_t ml_ode_backward_euler_sys (ml_sys_func_t f, double t0,
                                              const double *y0, int n, double dt,
                                              int steps, void *ctx,
                                              double *y_out) {
    double y [ML_ODE_MAXN], base [ML_ODE_MAXN];
    int i, s;
    if (ML_UNLIKELY (!f || !y0 || !y_out || n <= 0 || n > ML_ODE_MAXN))
        return ML_ERR_INVALID_ARG;
    for (i = 0; i < n; i++) y [i] = y0 [i];
    for (s = 0; s < steps; s++) {
        double tn = t0 + dt * (double) (s + 1);
        for (i = 0; i < n; i++) base [i] = y [i];
        if (oe_newton_step (f, tn, base, n, dt, ctx, y) != ML_SUCCESS)
            return ML_ERR_SINGULAR;
    }
    for (i = 0; i < n; i++) y_out [i] = y [i];
    return ML_SUCCESS;
}
ML_API ml_status_t ml_ode_trapezoidal_sys (ml_sys_func_t f, double t0,
                                           const double *y0, int n, double dt,
                                           int steps, void *ctx, double *y_out) {
    double y [ML_ODE_MAXN], base [ML_ODE_MAXN], fn [ML_ODE_MAXN];
    int i, s;
    if (ML_UNLIKELY (!f || !y0 || !y_out || n <= 0 || n > ML_ODE_MAXN))
        return ML_ERR_INVALID_ARG;
    for (i = 0; i < n; i++) y [i] = y0 [i];
    for (s = 0; s < steps; s++) {
        double tn = t0 + dt * (double) s, tn1 = tn + dt;
        f (tn, y, fn, n, ctx);
        for (i = 0; i < n; i++) base [i] = y [i] + 0.5 * dt * fn [i];
        if (oe_newton_step (f, tn1, base, n, 0.5 * dt, ctx, y) != ML_SUCCESS)
            return ML_ERR_SINGULAR;
    }
    for (i = 0; i < n; i++) y_out [i] = y [i];
    return ML_SUCCESS;
}
ML_API ml_status_t ml_ode_bdf2_sys (ml_sys_func_t f, double t0,
                                    const double *y0, int n, double dt,
                                    int steps, void *ctx, double *y_out) {
    double ym1 [ML_ODE_MAXN], y [ML_ODE_MAXN], base [ML_ODE_MAXN];
    int i, s;
    if (ML_UNLIKELY (!f || !y0 || !y_out || n <= 0 || n > ML_ODE_MAXN))
        return ML_ERR_INVALID_ARG;
    for (i = 0; i < n; i++) y [i] = y0 [i];
    if (steps <= 0) {
        for (i = 0; i < n; i++) y_out [i] = y [i];
        return ML_SUCCESS;
    }
    /* startup with backward Euler */
    for (i = 0; i < n; i++) {
        ym1 [i] = y [i];
        base [i] = y [i];
    }
    if (oe_newton_step (f, t0 + dt, base, n, dt, ctx, y) != ML_SUCCESS)
        return ML_ERR_SINGULAR;
    for (s = 1; s < steps; s++) {
        double tn1 = t0 + dt * (double) (s + 1);
        for (i = 0; i < n; i++)
            base [i] = (4.0 * y [i] - ym1 [i]) / 3.0;
        for (i = 0; i < n; i++) ym1 [i] = y [i];
        if (oe_newton_step (f, tn1, base, n, 2.0 / 3.0 * dt, ctx, y) != ML_SUCCESS)
            return ML_ERR_SINGULAR;
    }
    for (i = 0; i < n; i++) y_out [i] = y [i];
    return ML_SUCCESS;
}
static const double ab_coef [5][4] = {
    {0, 0, 0, 0},
    {1.0, 0, 0, 0},
    {1.5, -0.5, 0, 0},
    {23.0 / 12.0, -16.0 / 12.0, 5.0 / 12.0, 0},
    {55.0 / 24.0, -59.0 / 24.0, 37.0 / 24.0, -9.0 / 24.0}};
ML_API ml_status_t ml_ode_adams_bashforth_sys (ml_sys_func_t f, double t0,
                                               const double *y0, int n,
                                               int order, double dt, int steps,
                                               void *ctx, double *y_out) {
    double y [ML_ODE_MAXN], fh [4][ML_ODE_MAXN], k1 [ML_ODE_MAXN],
        k2 [ML_ODE_MAXN], k3 [ML_ODE_MAXN], k4 [ML_ODE_MAXN], tmp [ML_ODE_MAXN];
    int i, s, j, have = 1;
    if (ML_UNLIKELY (!f || !y0 || !y_out || n <= 0 || n > ML_ODE_MAXN || order < 1 || order > 4))
        return ML_ERR_INVALID_ARG;
    for (i = 0; i < n; i++) y [i] = y0 [i];
    f (t0, y, fh [0], n, ctx);
    for (s = 0; s < steps; s++) {
        double t = t0 + dt * (double) s;
        if (have < order) {
            double h2 = 0.5 * dt;
            f (t, y, k1, n, ctx);
            for (i = 0; i < n; i++) tmp [i] = y [i] + h2 * k1 [i];
            f (t + h2, tmp, k2, n, ctx);
            for (i = 0; i < n; i++) tmp [i] = y [i] + h2 * k2 [i];
            f (t + h2, tmp, k3, n, ctx);
            for (i = 0; i < n; i++) tmp [i] = y [i] + dt * k3 [i];
            f (t + dt, tmp, k4, n, ctx);
            for (i = 0; i < n; i++)
                y [i] += dt / 6.0 * (k1 [i] + 2.0 * k2 [i] + 2.0 * k3 [i] + k4 [i]);
        } else {
            for (i = 0; i < n; i++) {
                double acc = 0.0;
                for (j = 0; j < order; j++) acc += ab_coef [order][j] * fh [j][i];
                y [i] += dt * acc;
            }
        }
        for (j = 3; j > 0; j--)
            for (i = 0; i < n; i++) fh [j][i] = fh [j - 1][i];
        f (t + dt, y, fh [0], n, ctx);
        if (have < 4) have++;
    }
    for (i = 0; i < n; i++) y_out [i] = y [i];
    return ML_SUCCESS;
}
ML_API ml_status_t ml_ode_adams_moulton_sys (ml_sys_func_t f, double t0,
                                             const double *y0, int n, double dt,
                                             int steps, void *ctx,
                                             double *y_out) {
    double y [ML_ODE_MAXN], fn [ML_ODE_MAXN], fm1 [ML_ODE_MAXN],
        yp [ML_ODE_MAXN], fp [ML_ODE_MAXN];
    int i, s, it, have = 0;
    if (ML_UNLIKELY (!f || !y0 || !y_out || n <= 0 || n > ML_ODE_MAXN))
        return ML_ERR_INVALID_ARG;
    for (i = 0; i < n; i++) y [i] = y0 [i];
    for (s = 0; s < steps; s++) {
        double t = t0 + dt * (double) s;
        f (t, y, fn, n, ctx);
        if (!have) {
            /* Euler start */
            for (i = 0; i < n; i++) y [i] += dt * fn [i];
        } else {
            /* AB2 predictor */
            for (i = 0; i < n; i++)
                yp [i] = y [i] + dt * (1.5 * fn [i] - 0.5 * fm1 [i]);
            for (it = 0; it < 2; it++) {
                f (t + dt, yp, fp, n, ctx);
                for (i = 0; i < n; i++)
                    yp [i] = y [i] + 0.5 * dt * (fp [i] + fn [i]);
            }
            for (i = 0; i < n; i++) y [i] = yp [i];
        }
        for (i = 0; i < n; i++) fm1 [i] = fn [i];
        have = 1;
    }
    for (i = 0; i < n; i++) y_out [i] = y [i];
    return ML_SUCCESS;
}
/* ------------------------------------------------------------- BVP ------- */

ML_API ml_status_t ml_bvp_fd_linear (ml_func_t p, ml_func_t q, ml_func_t r,
                                     double x0, double x1, double ya, double yb,
                                     int n, double *x, double *y) {
    double a [512], b [512], c [512], d [512], h;
    int i;
    if (ML_UNLIKELY (!x || !y || n <= 0 || n > 510 || x1 <= x0))
        return ML_ERR_INVALID_ARG;
    h = (x1 - x0) / (double) (n + 1);
    for (i = 0; i <= n + 1; i++) x [i] = x0 + h * (double) i;
    y [0] = ya;
    y [n + 1] = yb;
    for (i = 1; i <= n; i++) {
        double xi = x [i], pi = p ? p (xi) : 0.0;
        double qi = q ? q (xi) : 0.0, ri = r ? r (xi) : 0.0;
        a [i - 1] = 1.0 / (h * h) - pi / (2.0 * h);
        b [i - 1] = -2.0 / (h * h) + qi;
        c [i - 1] = 1.0 / (h * h) + pi / (2.0 * h);
        d [i - 1] = ri;
    }
    d [0] -= a [0] * ya;
    d [n - 1] -= c [n - 1] * yb;
    /* Thomas algorithm */
    for (i = 1; i < n; i++) {
        double m = a [i] / b [i - 1];
        b [i] -= m * c [i - 1];
        d [i] -= m * d [i - 1];
    }
    y [n] = d [n - 1] / b [n - 1];
    for (i = n - 2; i >= 0; i--) y [i + 1] = (d [i] - c [i] * y [i + 2]) / b [i];
    return ML_SUCCESS;
}
/* -------------------------------------------------------- classics ------- */

ML_API void ml_lorenz_rhs (double t, const double *y, double *dy, int n,
                           void *ctx) {
    double sigma = 10.0, rho = 28.0, beta = 8.0 / 3.0;
    (void) t;
    if (n < 3) return;
    if (ctx) {
        const double *p = (const double *) ctx;
        sigma = p [0];
        rho = p [1];
        beta = p [2];
    }
    dy [0] = sigma * (y [1] - y [0]);
    dy [1] = y [0] * (rho - y [2]) - y [1];
    dy [2] = y [0] * y [1] - beta * y [2];
}
ML_API void ml_rossler_rhs (double t, const double *y, double *dy, int n,
                            void *ctx) {
    double a = 0.2, b = 0.2, c = 5.7;
    (void) t;
    if (n < 3) return;
    if (ctx) {
        const double *p = (const double *) ctx;
        a = p [0];
        b = p [1];
        c = p [2];
    }
    dy [0] = -y [1] - y [2];
    dy [1] = y [0] + a * y [1];
    dy [2] = b + y [2] * (y [0] - c);
}
