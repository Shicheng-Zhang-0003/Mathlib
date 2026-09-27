#include "ml_compiler.h"
#include "ml_optim_n.h"
#define ML_ON_MAX 32
static double ml_on_fdgrad(ml_vec_func_t f, const double *x, int n, int i) {
    double xs[32], xp[32], xm[32];
    for (int k = 0; k < n; k++) { xs[k] = x[k]; xp[k] = x[k]; xm[k] = x[k]; }
    (void)xs;
    double h = 1.4901161193847656e-08 * (1.0 + ml_fabs(x[i]));
    xp[i] = x[i] + h; xm[i] = x[i] - h;
    if (xp[i] == xm[i]) return ml_make_nan();
    double fp = f(xp, n), fm = f(xm, n);
    if (!ml_isfinite(fp) || !ml_isfinite(fm)) return ml_make_nan();
    return (fp - fm) / (2.0 * h);
}
ML_API ml_status_t ml_nelder_mead(ml_vec_func_t f, const double *x0, int n,
                                  double step, double tol, int max_iter,
                                  double *x_out, double *f_out) {
    if (!f || !x0 || !x_out || n <= 0 || n > 16) return ML_ERR_INVALID_ARG;
    if (!(step > 0.0) || !ml_isfinite(step)) step = 0.1;
    if (!(tol > 0.0) || !ml_isfinite(tol)) return ML_ERR_INVALID_ARG;
    if (max_iter <= 0) max_iter = 500;
    /* Stack-local simplex (never static): thread-safe, no cross-call state. */
    double S[17][16]; double F[17];
    for (int i = 0; i <= n; i++) {
        for (int j = 0; j < n; j++) S[i][j] = x0[j] + ((i > 0 && i - 1 == j) ? step * (1.0 + ml_fabs(x0[j])) : 0.0);
        F[i] = f(S[i], n);
        if (!ml_isfinite(F[i])) return ML_ERR_NAN_INPUT;
    }
    for (int it = 0; it < max_iter; it++) {
        int lo = 0, hi = 0;
        for (int i = 1; i <= n; i++) { if (F[i] < F[lo]) lo = i; if (F[i] > F[hi]) hi = i; }
        int hi2 = lo;
        for (int i = 0; i <= n; i++) { if (i != hi && F[i] > F[hi2]) hi2 = i; }
        double frange = ml_fabs(F[hi] - F[lo]);
        /* Simplex diameter max||S_i-S_lo||: require BOTH the function
         * range and the simplex footprint small. Range alone declares
         * victory on flat plateaus far from the minimizer. */
        double diam = 0.0;
        double xnorm = 0.0;
        {
            long double xn2 = 0.0L;
            for (int j = 0; j < n; j++) xn2 += (long double)S[lo][j] * (long double)S[lo][j];
            xnorm = (double)__builtin_sqrtl(xn2);
        }
        for (int i = 0; i <= n; i++) {
            if (i == lo) continue;
            long double s2 = 0.0L;
            for (int j = 0; j < n; j++) {
                long double dd = (long double)S[i][j] - (long double)S[lo][j];
                s2 += dd * dd;
            }
            {
                double dn = (double)__builtin_sqrtl(s2);
                if (dn > diam) diam = dn;
            }
        }
        if (frange <= tol * (1.0 + ml_fabs(F[lo])) &&
            diam <= tol * (1.0 + xnorm)) {
            for (int j = 0; j < n; j++) x_out[j] = S[lo][j];
            if (f_out) *f_out = F[lo];
            return ML_SUCCESS;
        }
        double C[16]; for (int j = 0; j < n; j++) { long double s = 0; for (int i = 0; i <= n; i++) if (i != hi) s += S[i][j]; C[j] = (double)(s / n); }
        double Xr[16]; for (int j = 0; j < n; j++) Xr[j] = C[j] + (C[j] - S[hi][j]);
        double fr = f(Xr, n); if (!ml_isfinite(fr)) return ML_ERR_SINGULAR;
        if (fr < F[lo]) {
            double Xe[16]; for (int j = 0; j < n; j++) Xe[j] = C[j] + 2.0 * (Xr[j] - C[j]);
            double fe = f(Xe, n); if (!ml_isfinite(fe)) return ML_ERR_SINGULAR;
            if (fe < fr) { for (int j = 0; j < n; j++) S[hi][j] = Xe[j]; F[hi] = fe; }
            else { for (int j = 0; j < n; j++) S[hi][j] = Xr[j]; F[hi] = fr; }
        } else if (fr < F[hi2]) { for (int j = 0; j < n; j++) S[hi][j] = Xr[j]; F[hi] = fr; }
        else {
            double Xc[16];
            if (fr < F[hi]) { for (int j = 0; j < n; j++) Xc[j] = C[j] + 0.5 * (Xr[j] - C[j]); }
            else { for (int j = 0; j < n; j++) Xc[j] = C[j] + 0.5 * (S[hi][j] - C[j]); }
            double fc = f(Xc, n); if (!ml_isfinite(fc)) return ML_ERR_SINGULAR;
            if (fc < F[hi]) { for (int j = 0; j < n; j++) S[hi][j] = Xc[j]; F[hi] = fc; }
            else { for (int i = 0; i <= n; i++) { if (i == lo) continue; for (int j = 0; j < n; j++) S[i][j] = S[lo][j] + 0.5 * (S[i][j] - S[lo][j]); F[i] = f(S[i], n); if (!ml_isfinite(F[i])) return ML_ERR_SINGULAR; } }
        }
    }
    return ML_ERR_SINGULAR;
}
ML_API ml_status_t ml_lbfgs_min(ml_vec_func_t f, const double *x0, int n,
                                double tol, int max_iter, double *x_out, double *f_out) {
    if (!f || !x0 || !x_out || n <= 0 || n > 32) return ML_ERR_INVALID_ARG;
    if (!(tol > 0.0) || !ml_isfinite(tol)) return ML_ERR_INVALID_ARG;
    if (max_iter <= 0) max_iter = 200;
    /* Stack locals (never static storage): thread-safe, no cross-call state. */
    double x[32], g[32], d[32], xn[32], gn[32];
    double sH[5][32], yH[5][32]; double rho[5];
    for (int j = 0; j < n; j++) x[j] = x0[j];
    double fx = f(x, n); if (!ml_isfinite(fx)) return ML_ERR_NAN_INPUT;
    for (int j = 0; j < n; j++) { g[j] = ml_on_fdgrad(f, x, n, j); if (!ml_isfinite(g[j])) return ML_ERR_SINGULAR; }
    /* Entry gradient norm: converge on ||g||<=tol*(1+||g0||), i.e.
     * gradient against gradient. Mixing in |fx| tests the wrong scale
     * (flat f with steep g, or steep f with tiny g, both mis-decide). */
    double g0n = 0.0;
    {
        long double s02 = 0.0L;
        for (int j = 0; j < n; j++) s02 += (long double)g[j] * (long double)g[j];
        g0n = (double)__builtin_sqrtl(s02);
    }
    int mem = 0, start = 0;
    for (int it = 0; it < max_iter; it++) {
        long double gn2 = 0; for (int j = 0; j < n; j++) gn2 += (long double)g[j] * g[j];
        if (!(gn2 > 0)) { for (int j = 0; j < n; j++) { x_out[j] = x[j]; }
            if (f_out) { *f_out = fx; }
            return ML_SUCCESS; }
        if (__builtin_sqrtl(gn2) <= (long double)tol * (1.0L + (long double)g0n)) {
            for (int j = 0; j < n; j++) { x_out[j] = x[j]; }
            if (f_out) { *f_out = fx; }
            return ML_SUCCESS;
        }
        /* two-loop recursion, double throughout: float casts truncated
         * the direction to ~1e-7 and stalled superlinear convergence. */
        double q[32], r[32], al[5];
        for (int j = 0; j < n; j++) q[j] = g[j];
        int m = mem < 5 ? mem : 5;
        for (int k = m - 1; k >= 0; k--) {
            int idx = (start - 1 - k + 10) % 5;
            al[k] = 0.0; for (int j = 0; j < n; j++) al[k] += sH[idx][j] * q[j];
            al[k] *= rho[idx];
            for (int j = 0; j < n; j++) q[j] -= al[k] * yH[idx][j];
        }
        double gamma = 1.0;
        if (mem > 0) {
            int idx = (start - 1 + 5) % 5;
            long double sy = 0, yy = 0;
            for (int j = 0; j < n; j++) { sy += (long double)sH[idx][j] * yH[idx][j]; yy += (long double)yH[idx][j] * yH[idx][j]; }
            if (yy > 0) gamma = (double)(sy / yy);
        }
        for (int j = 0; j < n; j++) r[j] = gamma * q[j];
        for (int k = 0; k < m; k++) {
            int idx = (start - m + k + 10) % 5;
            double beta = 0.0; for (int j = 0; j < n; j++) beta += yH[idx][j] * r[j];
            beta *= rho[idx];
            for (int j = 0; j < n; j++) r[j] += sH[idx][j] * (al[k] - beta);
        }
        for (int j = 0; j < n; j++) d[j] = -r[j];
        /* backtracking Armijo */
        double gd = 0; for (int j = 0; j < n; j++) gd += g[j] * d[j];
        if (!(gd < 0)) { for (int j = 0; j < n; j++) d[j] = -g[j]; gd = -(double)gn2; }
        double alpha = 1.0;
        for (int ls = 0; ls < 40; ls++) {
            for (int j = 0; j < n; j++) xn[j] = x[j] + alpha * d[j];
            double fn = f(xn, n); if (!ml_isfinite(fn)) { alpha *= 0.5; continue; }
            if (fn <= fx + 1e-4 * alpha * gd) { fx = fn; break; }
            alpha *= 0.5;
            if (alpha < 1e-16) return ML_ERR_SINGULAR;
        }
        for (int j = 0; j < n; j++) { gn[j] = ml_on_fdgrad(f, xn, n, j); if (!ml_isfinite(gn[j])) return ML_ERR_SINGULAR; }
        int idx = start % 5;
        for (int j = 0; j < n; j++) { sH[idx][j] = xn[j] - x[j]; yH[idx][j] = gn[j] - g[j]; }
        {
            long double sy = 0; for (int j = 0; j < n; j++) sy += (long double)sH[idx][j] * yH[idx][j];
            if (sy <= 0) { /* skip update */ }
            else { rho[idx] = (double)(1.0L / sy); start++; mem++; }
        }
        for (int j = 0; j < n; j++) { x[j] = xn[j]; g[j] = gn[j]; }
    }
    return ML_ERR_SINGULAR;
}
ML_API ml_status_t ml_adam_min(ml_vec_func_t f, const double *x0, int n,
                               double lr, double tol, int max_iter, double *x_out, double *f_out) {
    if (!f || !x0 || !x_out || n <= 0 || n > 32) return ML_ERR_INVALID_ARG;
    if (!(lr > 0.0) || !ml_isfinite(lr)) return ML_ERR_INVALID_ARG;
    if (!(tol > 0.0) || !ml_isfinite(tol)) return ML_ERR_INVALID_ARG;
    if (max_iter <= 0) max_iter = 2000;
    /* Stack locals (never static storage): thread-safe, no cross-call state. */
    double x[32], m[32], v[32], g[32];
    for (int j = 0; j < n; j++) { x[j] = x0[j]; m[j] = 0.0; v[j] = 0.0; }
    double fx = f(x, n); if (!ml_isfinite(fx)) return ML_ERR_NAN_INPUT;
    for (int t = 1; t <= max_iter; t++) {
        /* Full gradient at the CURRENT x first: evaluate-and-step per
         * coordinate (Gauss-Seidel) differentiates a moving point, biasing
         * the m/v moments and breaking Adam's semantics. */
        double gmax = 0.0;
        for (int j = 0; j < n; j++) {
            g[j] = ml_on_fdgrad(f, x, n, j);
            if (!ml_isfinite(g[j])) return ML_ERR_SINGULAR;
            {
                double ag = ml_fabs(g[j]);
                if (ag > gmax) gmax = ag;
            }
        }
        {
            double bc1 = 1.0 - __builtin_pow(0.9, t);
            double bc2 = 1.0 - __builtin_pow(0.999, t);
            for (int j = 0; j < n; j++) {
                m[j] = 0.9 * m[j] + 0.1 * g[j];
                v[j] = 0.999 * v[j] + 0.001 * g[j] * g[j];
                {
                    double mh = m[j] / bc1;
                    double vh = v[j] / bc2;
                    double dx = lr * mh / (__builtin_sqrt(vh) + 1e-8);
                    x[j] -= dx;
                }
            }
        }
        double fn = f(x, n); if (!ml_isfinite(fn)) return ML_ERR_SINGULAR;
        fx = fn;
        if (gmax <= tol) { for (int j = 0; j < n; j++) { x_out[j] = x[j]; }
            if (f_out) { *f_out = fx; }
            return ML_SUCCESS; }
    }
    return ML_ERR_SINGULAR;
}
