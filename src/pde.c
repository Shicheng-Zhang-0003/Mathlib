#include "ml_compiler.h"
#include "ml_pde.h"
ML_API ml_status_t ml_heat_explicit(const double *u0, double *u1, int n, double dx, double dt, double alpha) {
    if (!u0 || !u1 || n < 3 || n > 4096) return ML_ERR_INVALID_ARG;
    if (!(dx > 0.0) || !(dt > 0.0) || !(alpha > 0.0)) return ML_ERR_INVALID_ARG;
    if (!ml_isfinite(dx) || !ml_isfinite(dt) || !ml_isfinite(alpha)) return ML_ERR_NAN_INPUT;
    double r = alpha * dt / (dx * dx);
    if (!(r <= 0.5)) return ML_ERR_SINGULAR;
    for (int i = 0; i < n; i++) if (!ml_isfinite(u0[i])) return ML_ERR_NAN_INPUT;
    u1[0] = u0[0]; u1[n-1] = u0[n-1];
    for (int i = 1; i < n - 1; i++) {
        u1[i] = u0[i] + r * (u0[i+1] - 2.0 * u0[i] + u0[i-1]);
        if (!ml_isfinite(u1[i])) return ML_ERR_SINGULAR;
    }
    return ML_SUCCESS;
}
ML_API ml_status_t ml_heat_implicit(const double *u0, double *u1, int n, double dx, double dt, double alpha) {
    if (!u0 || !u1 || n < 3 || n > 1024) return ML_ERR_INVALID_ARG;
    if (!(dx > 0.0) || !(dt > 0.0) || !(alpha > 0.0)) return ML_ERR_INVALID_ARG;
    double r = alpha * dt / (dx * dx);
    if (!ml_isfinite(r)) return ML_ERR_NAN_INPUT;
    for (int i = 0; i < n; i++) if (!ml_isfinite(u0[i])) return ML_ERR_NAN_INPUT;
    /* (I + r*T) u1 = u0 with Dirichlet ends; Thomas. */
    double cp[1024], dp[1024];
    u1[0] = u0[0]; u1[n-1] = u0[n-1];
    if (n == 3) { u1[1] = (u0[1] + r * (u0[0] + u0[2])) / (1.0 + 2.0 * r); return ML_SUCCESS; }
    {
        double den = 1.0 + 2.0 * r;
        cp[1] = -r / den;
        dp[1] = (u0[1] + r * u0[0]) / den;
        for (int i = 2; i < n - 1; i++) {
            double rhs = u0[i] + (i == n - 2 ? r * u0[n-1] : 0.0);
            double d2 = 1.0 + 2.0 * r + r * cp[i-1];
            if (d2 == 0.0) return ML_ERR_SINGULAR;
            cp[i] = -r / d2;
            dp[i] = (rhs + r * dp[i-1]) / d2;
        }
        u1[n-2] = dp[n-2];
        for (int i = n - 3; i >= 1; i--) u1[i] = dp[i] - cp[i] * u1[i+1];
        return ML_SUCCESS;
    }
}
ML_API ml_status_t ml_poisson_1d(const double *f, double *u, int n, double dx) {
    if (!f || !u || n < 3 || n > 1024) return ML_ERR_INVALID_ARG;
    if (!(dx > 0.0) || !ml_isfinite(dx)) return ML_ERR_INVALID_ARG;
    for (int i = 0; i < n; i++) if (!ml_isfinite(f[i])) return ML_ERR_NAN_INPUT;
    /* -u''=f, Dirichlet 0 ends: -u_{i-1}+2u_i-u_{i+1}=dx^2 f_i; Thomas. */
    double cp[1024], dp[1024];
    u[0] = 0.0; u[n-1] = 0.0;
    if (n == 3) { u[1] = 0.5 * dx * dx * f[1]; return ML_SUCCESS; }
    {
        double h2 = dx * dx;
        cp[1] = -1.0 / 2.0;
        dp[1] = h2 * f[1] / 2.0;
        for (int i = 2; i < n - 1; i++) {
            double den = 2.0 + cp[i-1];
            if (den == 0.0) return ML_ERR_SINGULAR;
            cp[i] = -1.0 / den;
            dp[i] = (h2 * f[i] + dp[i-1]) / den;
        }
        u[n-2] = dp[n-2];
        for (int i = n - 3; i >= 1; i--) u[i] = dp[i] - cp[i] * u[i+1];
        return ML_SUCCESS;
    }
}
ML_API ml_status_t ml_wave_leapfrog(const double *u_prev, const double *u_cur, double *u_next, int n, double cfl) {
    if (!u_prev || !u_cur || !u_next || n < 3 || n > 4096) return ML_ERR_INVALID_ARG;
    if (!(cfl >= 0.0) || !(cfl <= 1.0) || !ml_isfinite(cfl)) return ML_ERR_SINGULAR;
    double c2 = cfl * cfl;
    u_next[0] = u_cur[0]; u_next[n-1] = u_cur[n-1];
    for (int i = 1; i < n - 1; i++) {
        if (!ml_isfinite(u_prev[i]) || !ml_isfinite(u_cur[i]) || !ml_isfinite(u_cur[i+1]) || !ml_isfinite(u_cur[i-1])) return ML_ERR_NAN_INPUT;
        u_next[i] = 2.0 * u_cur[i] - u_prev[i] + c2 * (u_cur[i+1] - 2.0 * u_cur[i] + u_cur[i-1]);
        if (!ml_isfinite(u_next[i])) return ML_ERR_SINGULAR;
    }
    return ML_SUCCESS;
}
ML_API ml_status_t ml_fem1d_assemble(double h, double *ke00, double *ke01, double *ke11, double *me00, double *me01, double *me11) {
    if (!ke00 || !ke01 || !ke11 || !me00 || !me01 || !me11) return ML_ERR_INVALID_ARG;
    if (!(h > 0.0) || !ml_isfinite(h)) return ML_ERR_INVALID_ARG;
    *ke00 = 1.0 / h; *ke01 = -1.0 / h; *ke11 = 1.0 / h;
    *me00 = h / 3.0; *me01 = h / 6.0; *me11 = h / 3.0;
    return ML_SUCCESS;
}
