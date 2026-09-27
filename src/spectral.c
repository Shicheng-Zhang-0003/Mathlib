#include "ml_compiler.h"
#include "ml_spectral.h"
ML_API ml_status_t ml_cg_solve(ml_tensor_view_t A, const double *b, double *x, int max_iter, double tol) {
    int n = A.rows;
    if (!A.data || !b || !x || n <= 0 || n > 256) return ML_ERR_INVALID_ARG;
    if (A.cols != n) return ML_ERR_INVALID_ARG;
    if (!(tol > 0.0) || !ml_isfinite(tol)) return ML_ERR_INVALID_ARG;
    if (max_iter <= 0) max_iter = 4 * n;
    /* Stack-local scratch (thread-safe; no static storage). */
    double r[256], p[256], Ap[256];
    for (int i = 0; i < n; i++) {
        if (!ml_isfinite(b[i])) return ML_ERR_NAN_INPUT;
        for (int j = 0; j < n; j++) if (!ml_isfinite(ML_TENSOR_AT(A,i,j))) return ML_ERR_NAN_INPUT;
        x[i] = 0.0; r[i] = b[i]; p[i] = b[i];
    }
    long double rs = 0; for (int i = 0; i < n; i++) rs += (long double)r[i]*r[i];
    if (!(rs > 0)) return ML_SUCCESS;
    for (int it = 0; it < max_iter; it++) {
        for (int i = 0; i < n; i++) { long double s = 0; for (int j = 0; j < n; j++) s += (long double)ML_TENSOR_AT(A,i,j)*p[j]; Ap[i] = (double)s; }
        long double pAp = 0; for (int i = 0; i < n; i++) pAp += (long double)p[i]*Ap[i];
        if (!(pAp > 0)) return ML_ERR_SINGULAR;
        double alpha = (double)(rs / pAp);
        long double rsn = 0;
        for (int i = 0; i < n; i++) { x[i] += alpha * p[i]; r[i] -= alpha * Ap[i]; rsn += (long double)r[i]*r[i]; }
        if (__builtin_sqrtl(rsn) <= (long double)tol) return ML_SUCCESS;
        double beta = (double)(rsn / rs);
        for (int i = 0; i < n; i++) p[i] = r[i] + beta * p[i];
        rs = rsn;
    }
    return ML_ERR_SINGULAR;
}
ML_API ml_status_t ml_gmres_solve(ml_tensor_view_t A, const double *b, double *x, int restart, int max_iter, double tol) {
    int n = A.rows;
    if (!A.data || !b || !x || n <= 0 || n > 64) return ML_ERR_INVALID_ARG;
    if (A.cols != n) return ML_ERR_INVALID_ARG;
    if (restart <= 0) restart = n < 20 ? n : 20;
    if (restart > 32) restart = 32;
    if (!(tol > 0.0) || !ml_isfinite(tol)) return ML_ERR_INVALID_ARG;
    if (max_iter <= 0) max_iter = 4 * n;
    for (int i = 0; i < n; i++) x[i] = 0.0;
    /* Stack-local scratch (thread-safe; no static storage). */
    double r[64], V[33][64], H[33][32], cs[32], sn[32], g[33], ys[32];
    for (int outer = 0; outer < max_iter; outer++) {
        for (int i = 0; i < n; i++) { long double s = (long double)b[i]; for (int j = 0; j < n; j++) s -= (long double)ML_TENSOR_AT(A,i,j)*x[j]; r[i] = (double)s; }
        long double beta = 0; for (int i = 0; i < n; i++) beta += (long double)r[i]*r[i];
        beta = __builtin_sqrtl(beta);
        if (beta <= (long double)tol) return ML_SUCCESS;
        for (int i = 0; i < n; i++) V[0][i] = (double)((long double)r[i]/beta);
        for (int i = 0; i <= restart; i++) { g[i] = 0; }
        g[0] = (double)beta;
        int j = 0;
        for (; j < restart; j++) {
            double w[64];
            for (int i = 0; i < n; i++) { long double s = 0; for (int k = 0; k < n; k++) s += (long double)ML_TENSOR_AT(A,i,k)*V[j][k]; w[i] = (double)s; }
            for (int i = 0; i <= j; i++) {
                long double h = 0; for (int k = 0; k < n; k++) h += (long double)w[k]*V[i][k];
                H[i][j] = (double)h;
                for (int k = 0; k < n; k++) w[k] -= H[i][j] * V[i][k];
            }
            long double nw = 0; for (int k = 0; k < n; k++) nw += (long double)w[k]*w[k];
            H[j+1][j] = (double)__builtin_sqrtl(nw);
            if (H[j+1][j] != 0) for (int k = 0; k < n; k++) V[j+1][k] = (double)((long double)w[k]/H[j+1][j]);
            /* Givens */
            for (int i = 0; i < j; i++) {
                double t1 = cs[i]*H[i][j] + sn[i]*H[i+1][j];
                double t2 = -sn[i]*H[i][j] + cs[i]*H[i+1][j];
                H[i][j] = t1; H[i+1][j] = t2;
            }
            double den = __builtin_sqrt(H[j][j]*H[j][j] + H[j+1][j]*H[j+1][j]);
            if (den == 0) break;
            cs[j] = H[j][j]/den; sn[j] = H[j+1][j]/den;
            H[j][j] = den; H[j+1][j] = 0;
            g[j+1] = -sn[j]*g[j]; g[j] *= cs[j];
            if (ml_fabs(g[j+1]) <= tol) { j++; break; }
        }
        for (int i = j - 1; i >= 0; i--) {
            ys[i] = g[i];
            for (int k = i + 1; k < j; k++) ys[i] -= H[i][k]*ys[k];
            if (H[i][i] == 0) return ML_ERR_SINGULAR;
            ys[i] /= H[i][i];
        }
        for (int i = 0; i < j; i++) for (int k = 0; k < n; k++) x[k] += ys[i]*V[i][k];
        if (ml_fabs(g[j]) <= tol) return ML_SUCCESS;
    }
    return ML_ERR_SINGULAR;
}
ML_API ml_status_t ml_power_iter(ml_tensor_view_t A, double *lambda, double *vec, int max_iter, double tol) {
    int n = A.rows;
    if (!A.data || !lambda || n <= 0 || n > 128) return ML_ERR_INVALID_ARG;
    if (A.cols != n) return ML_ERR_INVALID_ARG;
    if (max_iter <= 0) max_iter = 1000;
    if (!(tol > 0.0) || !ml_isfinite(tol)) tol = 1e-10;
    /* Stack-local scratch (thread-safe; no static storage). */
    double v[128], w[128];
    for (int i = 0; i < n; i++) v[i] = 1.0 / __builtin_sqrt((double)n);
    double lam = 0;
    for (int it = 0; it < max_iter; it++) {
        for (int i = 0; i < n; i++) { long double s = 0; for (int j = 0; j < n; j++) s += (long double)ML_TENSOR_AT(A,i,j)*v[j]; w[i] = (double)s; }
        long double nw = 0; for (int i = 0; i < n; i++) nw += (long double)w[i]*w[i];
        nw = __builtin_sqrtl(nw);
        if (!(nw > 0)) return ML_ERR_SINGULAR;
        for (int i = 0; i < n; i++) w[i] = (double)((long double)w[i]/nw);
        long double rq = 0;
        for (int i = 0; i < n; i++) { long double s = 0; for (int j = 0; j < n; j++) s += (long double)ML_TENSOR_AT(A,i,j)*w[j]; rq += (long double)w[i]*s; }
        if (ml_fabs((double)rq - lam) <= tol * (1.0 + ml_fabs((double)rq))) {
            *lambda = (double)rq;
            if (vec) for (int i = 0; i < n; i++) vec[i] = w[i];
            return ML_SUCCESS;
        }
        lam = (double)rq;
        for (int i = 0; i < n; i++) v[i] = w[i];
    }
    return ML_ERR_SINGULAR;
}
ML_API ml_status_t ml_svd_jacobi(ml_tensor_view_t A, double *svals, ml_tensor_view_t U, ml_tensor_view_t Vt, int max_sweeps) {
    int m = A.rows, n = A.cols;
    if (!A.data || !svals || m <= 0 || n <= 0 || m > 32 || n > 32) return ML_ERR_INVALID_ARG;
    if (max_sweeps <= 0) max_sweeps = 20;
    /* Stack-local scratch (thread-safe; no static storage). */
    double B[32][32], UU[32][32], VV[32][32];
    for (int i = 0; i < m; i++) for (int j = 0; j < n; j++) { double v = ML_TENSOR_AT(A,i,j); if (!ml_isfinite(v)) return ML_ERR_NAN_INPUT; B[i][j] = v; }
    for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) UU[i][j] = (i == j);
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) VV[i][j] = (i == j);
    int p = m < n ? m : n;
    for (int sw = 0; sw < max_sweeps; sw++) {
        long double off = 0;
        for (int i = 0; i < p; i++) for (int j = i + 1; j < n && j < 32; j++) {
            long double d = 0; for (int k = 0; k < m; k++) d += (long double)B[k][i]*B[k][j];
            off += d*d;
        }
        if (!(off > 1e-30L)) break;
        for (int i = 0; i < n - 1; i++) for (int j = i + 1; j < n; j++) {
            long double aii = 0, ajj = 0, aij = 0;
            for (int k = 0; k < m; k++) { aii += (long double)B[k][i]*B[k][i]; ajj += (long double)B[k][j]*B[k][j]; aij += (long double)B[k][i]*B[k][j]; }
            if ((double)aij == 0.0) continue;
            double tau = (double)((ajj - aii) / (2.0L * aij));
            double t = ((tau >= 0) ? 1.0 : -1.0) / (ml_fabs(tau) + __builtin_sqrt(tau*tau + 1.0));
            double c = 1.0 / __builtin_sqrt(t*t + 1.0), s = t * c;
            for (int k = 0; k < m; k++) { double bi = B[k][i], bj = B[k][j]; B[k][i] = c*bi - s*bj; B[k][j] = s*bi + c*bj; }
            for (int k = 0; k < n; k++) { double bi = VV[k][i], bj = VV[k][j]; VV[k][i] = c*bi - s*bj; VV[k][j] = s*bi + c*bj; }
        }
    }
    for (int j = 0; j < n; j++) {
        long double s = 0; for (int k = 0; k < m; k++) s += (long double)B[k][j]*B[k][j];
        svals[j] = (double)__builtin_sqrtl(s);
        if (svals[j] != 0 && U.data) for (int k = 0; k < m; k++) ML_TENSOR_AT(U,k,j) = (double)((long double)B[k][j]/svals[j]);
    }
    if (Vt.data) {
        if (Vt.rows != n || Vt.cols != n) return ML_ERR_INVALID_ARG;
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) ML_TENSOR_AT(Vt,i,j) = VV[j][i];
    }
    /* Sort singular values descending, permuting U columns and Vt rows
     * consistently via the same pairwise swaps (permutation tracking).
     * Previously only svals were swapped, leaving U/Vt mismatched. */
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (svals[j] > svals[i]) {
                double t = svals[i]; svals[i] = svals[j]; svals[j] = t;
                if (U.data) {
                    for (int k = 0; k < m; k++) {
                        double tu = ML_TENSOR_AT(U, k, i);
                        ML_TENSOR_AT(U, k, i) = ML_TENSOR_AT(U, k, j);
                        ML_TENSOR_AT(U, k, j) = tu;
                    }
                }
                if (Vt.data) {
                    for (int k = 0; k < n; k++) {
                        double tv = ML_TENSOR_AT(Vt, i, k);
                        ML_TENSOR_AT(Vt, i, k) = ML_TENSOR_AT(Vt, j, k);
                        ML_TENSOR_AT(Vt, j, k) = tv;
                    }
                }
            }
        }
    }
    (void)UU;
    return ML_SUCCESS;
}
ML_API ml_status_t ml_qr_iter_eig(ml_tensor_view_t A, double *evals_re, double *evals_im, int max_iter, double tol) {
    int n = A.rows;
    if (!A.data || !evals_re || !evals_im || n <= 0 || n > 16) return ML_ERR_INVALID_ARG;
    if (A.cols != n) return ML_ERR_INVALID_ARG;
    if (!(tol > 0.0) || !ml_isfinite(tol)) tol = 1e-10;
    if (max_iter <= 0) max_iter = 500;
    /* Stack-local scratch (thread-safe; no static storage). */
    double H[16][16];
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) { double v = ML_TENSOR_AT(A,i,j); if (!ml_isfinite(v)) return ML_ERR_NAN_INPUT; H[i][j] = v; }
    for (int it = 0; it < max_iter; it++) {
        long double off = 0;
        for (int i = 1; i < n; i++) for (int j = 0; j < i && j + 1 >= i; j++) { if (j + 1 == i || j + 1 < i) {} }
        for (int i = 1; i < n; i++) { long double v = H[i][i-1]; off += v*v; }
        if (__builtin_sqrtl(off) <= (long double)tol) break;
        /* Wilkinson shift from bottom 2x2 */
        double mu = H[n-1][n-1];
        if (n >= 2) {
            double a = H[n-2][n-2], b = H[n-2][n-1], c = H[n-1][n-2], d = H[n-1][n-1];
            double tr = a + d, det = a*d - b*c, disc = tr*tr - 4*det;
            if (disc >= 0) { double s = __builtin_sqrt(disc); double l1 = (tr+s)/2, l2 = (tr-s)/2; mu = (ml_fabs(l1-d) < ml_fabs(l2-d)) ? l1 : l2; }
        }
        for (int i = 0; i < n; i++) H[i][i] -= mu;
        /* explicit QR via Givens */
        for (int i = 0; i < n - 1; i++) {
            double x = H[i][i], y = H[i+1][i];
            if (y == 0) continue;
            double r = __builtin_sqrt(x*x + y*y);
            double c = x / r, s = -y / r;
            for (int k = i; k < n; k++) { double a1 = H[i][k], a2 = H[i+1][k]; H[i][k] = c*a1 - s*a2; H[i+1][k] = s*a1 + c*a2; }
            for (int k = 0; k < n; k++) { double a1 = H[k][i], a2 = H[k][i+1]; H[k][i] = c*a1 - s*a2; H[k][i+1] = s*a1 + c*a2; }
        }
        for (int i = 0; i < n; i++) H[i][i] += mu;
        if (it == max_iter - 1) return ML_ERR_SINGULAR;
    }
    for (int i = 0; i < n; i++) {
        if (i + 1 < n && ml_fabs(H[i+1][i]) > tol) {
            double a = H[i][i], b = H[i][i+1], c = H[i+1][i], d = H[i+1][i+1];
            double tr = a + d, det = a*d - b*c, disc = tr*tr - 4*det;
            if (disc < 0) { double s = __builtin_sqrt(-disc)/2; evals_re[i] = tr/2; evals_im[i] = s; evals_re[i+1] = tr/2; evals_im[i+1] = -s; i++; }
            else { evals_re[i] = a; evals_im[i] = 0; }
        } else { evals_re[i] = H[i][i]; evals_im[i] = 0; }
    }
    return ML_SUCCESS;
}
