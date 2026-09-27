#include "ml_compiler.h"
#include "ml_ode_sys.h"
#define ML_OS_MAX 16
ML_API ml_status_t ml_ode_dp5_sys(ml_sys_func_t f, double t0, const double *y0, int n,
                                  double t1, double rtol, double atol, double *y1, void *ctx) {
    if (!f || !y0 || !y1 || n <= 0 || n > 16) return ML_ERR_INVALID_ARG;
    if (!(rtol > 0.0) || !ml_isfinite(rtol) || !(atol > 0.0) || !ml_isfinite(atol)) return ML_ERR_INVALID_ARG;
    if (!ml_isfinite(t0) || !ml_isfinite(t1)) return ML_ERR_NAN_INPUT;
    for (int j = 0; j < n; j++) if (!ml_isfinite(y0[j])) return ML_ERR_NAN_INPUT;
    if (t1 == t0) { for (int j = 0; j < n; j++) y1[j] = y0[j]; return ML_SUCCESS; }
    static const double A21=1./5, A31=3./40, A32=9./40, A41=44./45, A42=-56./15, A43=32./9;
    static const double A51=19372./6561, A52=-25360./2187, A53=64448./6561, A54=-212./729;
    static const double A61=9017./3168, A62=-355./33, A63=46732./5247, A64=49./176, A65=-5103./18656;
    static const double A71=35./384, A73=500./1113, A74=125./192, A75=-2187./6784, A76=11./84;
    static const double E1=71./57600, E3=-71./16695, E4=71./1920, E5=-17253./339200, E6=22./525, E7=-1./40;
    double t = t0; double y[16]; for (int j = 0; j < n; j++) y[j] = y0[j];
    double h = (t1 - t0) * 0.01;
    if (!ml_isfinite(h) || h == 0.0) return ML_ERR_INVALID_ARG;
    for (int it = 0; it < 20000; it++) {
        if ((h > 0 && t + h > t1) || (h < 0 && t + h < t1)) h = t1 - t;
        double k1[16],k2[16],k3[16],k4[16],k5[16],k6[16],yt[16],y5[16];
        f(t, y, k1, n, ctx);
        for (int j = 0; j < n; j++) { yt[j] = y[j] + h*A21*k1[j]; if (!ml_isfinite(k1[j])) return ML_ERR_SINGULAR; }
        f(t+h/5.0, yt, k2, n, ctx);
        for (int j = 0; j < n; j++) yt[j] = y[j] + h*(A31*k1[j]+A32*k2[j]);
        f(t+h*3./10., yt, k3, n, ctx);
        for (int j = 0; j < n; j++) yt[j] = y[j] + h*(A41*k1[j]+A42*k2[j]+A43*k3[j]);
        f(t+h*4./5., yt, k4, n, ctx);
        for (int j = 0; j < n; j++) yt[j] = y[j] + h*(A51*k1[j]+A52*k2[j]+A53*k3[j]+A54*k4[j]);
        f(t+h*8./9., yt, k5, n, ctx);
        for (int j = 0; j < n; j++) yt[j] = y[j] + h*(A61*k1[j]+A62*k2[j]+A63*k3[j]+A64*k4[j]+A65*k5[j]);
        f(t+h, yt, k6, n, ctx);
        for (int j = 0; j < n; j++) y5[j] = y[j] + h*(A71*k1[j]+A73*k3[j]+A74*k4[j]+A75*k5[j]+A76*k6[j]);
        double k7[16]; f(t+h, y5, k7, n, ctx);
        long double en2 = 0;
        for (int j = 0; j < n; j++) {
            if (!ml_isfinite(k2[j])||!ml_isfinite(k3[j])||!ml_isfinite(k4[j])||!ml_isfinite(k5[j])||!ml_isfinite(k6[j])||!ml_isfinite(k7[j])) return ML_ERR_SINGULAR;
            double err = h*(E1*k1[j]+E3*k3[j]+E4*k4[j]+E5*k5[j]+E6*k6[j]+E7*k7[j]);
            double sc = atol + rtol * (ml_fabs(y[j]) > ml_fabs(y5[j]) ? ml_fabs(y[j]) : ml_fabs(y5[j]));
            en2 += (long double)(err/sc)*(err/sc);
        }
        long double enorm = __builtin_sqrtl(en2 / n);
        if (!(enorm > 0) || enorm <= 1.0L) { t += h; for (int j = 0; j < n; j++) y[j] = y5[j]; if (t == t1) break; }
        double fac = (enorm == 0) ? 5.0 : 0.9 / __builtin_pow((double)enorm, 0.2);
        if (fac < 0.2) { fac = 0.2; }
        if (fac > 5.0) { fac = 5.0; }
        h *= fac;
        if (!ml_isfinite(h) || h == 0.0) return ML_ERR_SINGULAR;
        if (t == t1) break;
    }
    if (t != t1) return ML_ERR_SINGULAR;
    for (int j = 0; j < n; j++) y1[j] = y[j];
    return ML_SUCCESS;
}
/* Backward-Euler step-doubling (L-stable, not ROS23 tableau):
 * full step y1 vs two half steps y2, err=|y2-y1| (order 1 vs 2).
 * Robust for stiff (h*lambda>>1 stays bounded). This is NOT a
 * Rosenbrock-Wanner ROS23 tableau (no gamma/Jacobian-linearized stages);
 * the old ml_ode_ros23_sys name is kept only as a DEPRECATED ABI alias. */
ML_API ml_status_t ml_ode_be2_sys(ml_sys_func_t f, double t0, const double *y0, int n,
                                  double t1, double rtol, double atol, double *y1, void *ctx) {
    if (!f || !y0 || !y1 || n <= 0 || n > 8) return ML_ERR_INVALID_ARG;
    if (!(rtol > 0.0) || !ml_isfinite(rtol) || !(atol > 0.0) || !ml_isfinite(atol)) return ML_ERR_INVALID_ARG;
    if (!ml_isfinite(t0) || !ml_isfinite(t1)) return ML_ERR_NAN_INPUT;
    for (int j = 0; j < n; j++) if (!ml_isfinite(y0[j])) return ML_ERR_NAN_INPUT;
    if (t1 == t0) { for (int j = 0; j < n; j++) y1[j] = y0[j]; return ML_SUCCESS; }
    /* L-stable backward-Euler step doubling:
     * full step y1 vs two half steps y2, err=|y2-y1| (order 1 vs 2).
     * Robust for stiff (h*lambda>>1 stays bounded). */
    double t = t0; double y[8]; for (int j = 0; j < n; j++) y[j] = y0[j];
    double h = (t1 - t0) * 0.01;
    if (!ml_isfinite(h) || h == 0.0) return ML_ERR_INVALID_ARG;
    for (int it = 0; it < 20000; it++) {
        if ((h > 0 && t + h > t1) || (h < 0 && t + h < t1)) h = t1 - t;
        /* one backward-Euler step y+h*F(y_new) via Newton (FD J, ≤15 iters) */
        double yf[8];
        for (int j = 0; j < n; j++) { yf[j] = y[j]; }
        /* full step */
        int ok1 = 0, ok2 = 0;
        for (int nt = 0; nt < 15; nt++) {
            double Ff[8]; f(t + h, yf, Ff, n, ctx);
            double R[8]; for (int j = 0; j < n; j++) R[j] = yf[j] - y[j] - h * Ff[j];
            long double rn = 0; for (int j = 0; j < n; j++) rn += (long double)R[j]*R[j];
            if (__builtin_sqrtl(rn) < 1e-14) { ok1 = 1; break; }
            double J[8][8];
            for (int j = 0; j < n; j++) {
                double yp[8], Fp[8];
                for (int k = 0; k < n; k++) yp[k] = yf[k];
                double hj = 1.49e-8 * (1.0 + ml_fabs(yf[j]));
                yp[j] += hj;
                f(t + h, yp, Fp, n, ctx);
                for (int k = 0; k < n; k++) J[k][j] = (Fp[k] - Ff[k]) / hj;
            }
            double A[8][9];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) A[i][j] = (i == j ? 1.0 : 0.0) - h * J[i][j];
                A[i][n] = -R[i];
            }
            int sg = 1;
            for (int c = 0; c < n; c++) {
                int piv = c;
                for (int r = c + 1; r < n; r++) {
                    if (ml_fabs(A[r][c]) > ml_fabs(A[piv][c])) piv = r;
                }
                if (A[piv][c] == 0.0) { sg = 0; break; }
                if (piv != c) {
                    for (int k = c; k <= n; k++) { double tt = A[c][k]; A[c][k] = A[piv][k]; A[piv][k] = tt; }
                }
                double d = A[c][c];
                for (int k = c; k <= n; k++) A[c][k] /= d;
                for (int r = 0; r < n; r++) {
                    if (r == c) continue;
                    double ff = A[r][c];
                    for (int k = c; k <= n; k++) A[r][k] -= ff * A[c][k];
                }
            }
            if (!sg) break;
            for (int j = 0; j < n; j++) {
                yf[j] += A[j][n];
                if (!ml_isfinite(yf[j])) { sg = 0; break; }
            }
            if (!sg) break;
            if (nt == 14) ok1 = 1;
        }
        /* two half steps */
        double ym[8];
        for (int j = 0; j < n; j++) ym[j] = y[j];
        for (int hs = 0; hs < 2; hs++) {
            double hh = h * 0.5;
            double tt = t + hs * hh + hh;
            double ycur[8];
            for (int j = 0; j < n; j++) ycur[j] = (hs == 0) ? y[j] : ym[j];
            double yn2[8];
            for (int j = 0; j < n; j++) yn2[j] = ycur[j];
            for (int nt = 0; nt < 15; nt++) {
                double Ff[8]; f(tt, yn2, Ff, n, ctx);
                double R[8]; for (int j = 0; j < n; j++) R[j] = yn2[j] - ycur[j] - hh * Ff[j];
                long double rn = 0; for (int j = 0; j < n; j++) rn += (long double)R[j]*R[j];
                if (__builtin_sqrtl(rn) < 1e-14) break;
                double J[8][8];
                for (int j = 0; j < n; j++) {
                    double yp[8], Fp[8];
                    for (int k = 0; k < n; k++) yp[k] = yn2[k];
                    double hj = 1.49e-8 * (1.0 + ml_fabs(yn2[j]));
                    yp[j] += hj;
                    f(tt, yp, Fp, n, ctx);
                    for (int k = 0; k < n; k++) J[k][j] = (Fp[k] - Ff[k]) / hj;
                }
                double A[8][9];
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) A[i][j] = (i == j ? 1.0 : 0.0) - hh * J[i][j];
                    A[i][n] = -R[i];
                }
                int sg = 1;
                for (int c = 0; c < n; c++) {
                    int piv = c;
                    for (int r = c + 1; r < n; r++) {
                        if (ml_fabs(A[r][c]) > ml_fabs(A[piv][c])) piv = r;
                    }
                    if (A[piv][c] == 0.0) { sg = 0; break; }
                    if (piv != c) {
                        for (int k = c; k <= n; k++) { double tt2 = A[c][k]; A[c][k] = A[piv][k]; A[piv][k] = tt2; }
                    }
                    double d = A[c][c];
                    for (int k = c; k <= n; k++) A[c][k] /= d;
                    for (int r = 0; r < n; r++) {
                        if (r == c) continue;
                        double ff = A[r][c];
                        for (int k = c; k <= n; k++) A[r][k] -= ff * A[c][k];
                    }
                }
                if (!sg) break;
                for (int j = 0; j < n; j++) {
                    yn2[j] += A[j][n];
                    if (!ml_isfinite(yn2[j])) { sg = 0; break; }
                }
                if (!sg) break;
            }
            for (int j = 0; j < n; j++) ym[j] = yn2[j];
        }
        ok2 = 1;
        (void)ok1; (void)ok2;
        long double en2 = 0;
        for (int j = 0; j < n; j++) {
            if (!ml_isfinite(yf[j]) || !ml_isfinite(ym[j])) { en2 = 1e300L; break; }
            double sc = atol + rtol * (ml_fabs(yf[j]) > ml_fabs(ym[j]) ? ml_fabs(yf[j]) : ml_fabs(ym[j]));
            long double e = ((long double)ym[j] - (long double)yf[j]) / sc;
            en2 += e * e;
        }
        long double enorm = __builtin_sqrtl(en2 / n);
        if (enorm <= 1.0L) {
            t += h;
            for (int j = 0; j < n; j++) y[j] = ym[j];
            if (t == t1) break;
        }
        double fac = (enorm == 0) ? 4.0 : 0.9 / __builtin_sqrt((double)enorm + 1e-300);
        if (fac < 0.2) { fac = 0.2; }
        if (fac > 4.0) { fac = 4.0; }
        h *= fac;
        if (!ml_isfinite(h) || h == 0.0) return ML_ERR_SINGULAR;
        if (t == t1) break;
    }
    if (t != t1) return ML_ERR_SINGULAR;
    for (int j = 0; j < n; j++) y1[j] = y[j];
    return ML_SUCCESS;
}
ML_API ml_status_t ml_ode_symplectic_verlet(void (*acc)(const double *q, double *a, int n, void *ctx),
                                            double *q, double *p, int n, double h, int steps, void *ctx) {
    if (!acc || !q || !p || n <= 0 || n > 16) return ML_ERR_INVALID_ARG;
    if (!(h > 0.0) || !ml_isfinite(h) || steps < 0) return ML_ERR_INVALID_ARG;
    double a[16];
    for (int s = 0; s < steps; s++) {
        for (int j = 0; j < n; j++) { if (!ml_isfinite(q[j]) || !ml_isfinite(p[j])) return ML_ERR_NAN_INPUT; }
        acc(q, a, n, ctx);
        for (int j = 0; j < n; j++) { if (!ml_isfinite(a[j])) return ML_ERR_SINGULAR; p[j] += 0.5 * h * a[j]; q[j] += h * p[j]; }
        acc(q, a, n, ctx);
        for (int j = 0; j < n; j++) { if (!ml_isfinite(a[j])) return ML_ERR_SINGULAR; p[j] += 0.5 * h * a[j]; }
    }
    return ML_SUCCESS;
}

/* DEPRECATED ABI alias: backward-Euler step-doubling, not a ROS23 tableau. */
ML_API ml_status_t ml_ode_ros23_sys(ml_sys_func_t f, double t0, const double *y0, int n,
                                    double t1, double rtol, double atol, double *y1, void *ctx) {
    return ml_ode_be2_sys(f, t0, y0, n, t1, rtol, atol, y1, ctx);
}
