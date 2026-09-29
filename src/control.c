#include "ml_compiler.h"
#include "ml_control.h"
ML_API ml_status_t ml_lqr_gain_2x2(double a00, double a01, double a10, double a11,
                                   double b0, double b1, double q, double r,
                                   double *k0, double *k1) {
    if (!k0 || !k1) return ML_ERR_INVALID_ARG;
    if (!ml_isfinite(a00)||!ml_isfinite(a01)||!ml_isfinite(a10)||!ml_isfinite(a11)) return ML_ERR_NAN_INPUT;
    if (!ml_isfinite(b0)||!ml_isfinite(b1)||!(q>0.0)||!(r>0.0)||!ml_isfinite(q)||!ml_isfinite(r)) return ML_ERR_INVALID_ARG;
    /* Solve continuous ARE A'P+PA-PBR^{-1}B'P+Q=0 via Newton-Kleinman
     * from P0=Q*I; 2x2 closed Lyapunov solves per step. */
    double p00 = q, p01 = 0, p11 = q;
    /* Newton-Kleinman needs A - G*P Hurwitz at every step, and P0 = Q*I
     * does not guarantee that: for the double integrator it gives
     * det(A_c) = 0, so the Lyapunov solve is singular and the iteration
     * aborts.  Nudge the initial P until the closed loop is stable. */
    {
        double g00 = b0 * b0 / r, g01 = b0 * b1 / r;
        double g10 = b1 * b0 / r, g11 = b1 * b1 / r;
        for (int attempt = 0; attempt < 500; attempt++) {
            double a00c = a00 - (g00 * p00 + g01 * p01);
            double a01c = a01 - (g00 * p01 + g01 * p11);
            double a10c = a10 - (g10 * p00 + g11 * p01);
            double a11c = a11 - (g10 * p01 + g11 * p11);
            double tr = a00c + a11c;
            double det = a00c * a11c - a01c * a10c;
            if (tr < 0.0 && det > 0.0) break;
            if (det <= 0.0) p01 += 0.25 * (1.0 + ml_fabs(p01));
            if (tr >= 0.0) { p00 += 0.5; p11 += 0.5; }
        }
    }
    for (int it = 0; it < 100; it++) {
        double s0 = (b0*b0*p00 + 2*b0*b1*p01 + b1*b1*p11) / r;
        double a00c = a00 - (b0*b0*p00 + b0*b1*p01) / r;
        double a01c = a01 - (b0*b0*p01 + b0*b1*p11) / r;
        double a10c = a10 - (b1*b0*p00 + b1*b1*p01) / r;
        double a11c = a11 - (b1*b0*p01 + b1*b1*p11) / r;
        (void)s0;
        /* Solve A_c'P+P A_c = -(Q+S) with S=PBR^{-1}B'P */
        double S00 = (b0*p00+b1*p01)*(b0*p00+b1*p01)/r;
        double S01 = (b0*p00+b1*p01)*(b0*p01+b1*p11)/r;
        double S11 = (b0*p01+b1*p11)*(b0*p01+b1*p11)/r;
        double q00 = -(q + S00), q01 = -S01, q11 = -(q + S11);
        /* Kronecker 3x3 for symmetric P */
        double M[3][4];
        M[0][0]=2*a00c; M[0][1]=2*a10c; M[0][2]=0; M[0][3]=q00;
        M[1][0]=a01c; M[1][1]=a00c+a11c; M[1][2]=a10c; M[1][3]=q01;
        M[2][0]=0; M[2][1]=2*a01c; M[2][2]=2*a11c; M[2][3]=q11;
        for (int c = 0; c < 3; c++) {
            int piv = c; for (int rr = c+1; rr < 3; rr++) if (ml_fabs(M[rr][c]) > ml_fabs(M[piv][c])) piv = rr;
            if (M[piv][c] == 0.0) return ML_ERR_SINGULAR;
            if (piv != c) for (int k = c; k < 4; k++) { double t=M[c][k]; M[c][k]=M[piv][k]; M[piv][k]=t; }
            double d = M[c][c]; for (int k = c; k < 4; k++) M[c][k] /= d;
            for (int rr = 0; rr < 3; rr++) { if (rr==c) continue; double ff=M[rr][c]; for (int k=c;k<4;k++) M[rr][k]-=ff*M[c][k]; }
        }
        double n00=M[0][3], n01=M[1][3], n11=M[2][3];
        if (!ml_isfinite(n00)||!ml_isfinite(n01)||!ml_isfinite(n11)) return ML_ERR_SINGULAR;
        double e = ml_fabs(n00-p00)+ml_fabs(n01-p01)+ml_fabs(n11-p11);
        p00=n00; p01=n01; p11=n11;
        if (e < 1e-12*(1.0+ml_fabs(p00)+ml_fabs(p11))) break;
        /* No dedicated non-convergence code in ml_status_t; the codebase
         * convention (spectral, ode_sys) reports iteration exhaustion as
         * ML_ERR_SINGULAR. Never return success without converging. */
        if (it == 99) return ML_ERR_SINGULAR;
    }
    *k0 = (b0*p00 + b1*p01)/r;
    *k1 = (b0*p01 + b1*p11)/r;
    return ML_SUCCESS;
}
ML_API ml_status_t ml_kalman_1d(double x0, double p0, const double *zs, int n,
                                double F, double H, double Q, double R,
                                double *x_out, double *p_out) {
    if (!zs || n <= 0 || n > 100000) return ML_ERR_INVALID_ARG;
    if (!ml_isfinite(x0)||!ml_isfinite(p0)||!(p0>=0.0)) return ML_ERR_NAN_INPUT;
    if (!ml_isfinite(F)||!ml_isfinite(H)||!(Q>=0.0)||!(R>0.0)) return ML_ERR_INVALID_ARG;
    double x = x0, p = p0;
    for (int i = 0; i < n; i++) {
        if (!ml_isfinite(zs[i])) return ML_ERR_NAN_INPUT;
        double xp = F * x, pp = F * F * p + Q;
        double y = zs[i] - H * xp;
        double S = H * H * pp + R;
        if (!(S > 0.0)) return ML_ERR_SINGULAR;
        double K = pp * H / S;
        x = xp + K * y;
        p = (1.0 - K * H) * pp;
        if (!ml_isfinite(x) || !ml_isfinite(p)) return ML_ERR_SINGULAR;
    }
    if (x_out) *x_out = x;
    if (p_out) *p_out = p;
    return ML_SUCCESS;
}
ML_API double ml_lyapunov_2x2_trace(double a00, double a01, double a10, double a11) {
    if (!ml_isfinite(a00)||!ml_isfinite(a01)||!ml_isfinite(a10)||!ml_isfinite(a11)) return ml_make_nan();
    /* tr(A) < 0 and det(A) > 0 iff 2x2 Hurwitz stable. Return margin. */
    (void)a00; (void)a01; (void)a10; (void)a11;
    double tr = a00 + a11, det = a00*a11 - a01*a10;
    if (tr < 0.0 && det > 0.0) return -tr;
    return ml_make_nan();
}
