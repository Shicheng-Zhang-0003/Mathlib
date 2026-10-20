#include "ml_compiler.h"
#include "ml_spectral.h"
ML_API ml_status_t ml_cg_solve (ml_tensor_view_t A, const double *b, double *x, int max_iter, double tol) {
    int n = A.rows;
    if (!A.data || !b || !x || n <= 0 || n > 256) return ML_ERR_INVALID_ARG;
    if (A.cols != n) return ML_ERR_INVALID_ARG;
    if (!(tol > 0.0) || !ml_isfinite (tol)) return ML_ERR_INVALID_ARG;
    if (max_iter <= 0) max_iter = 4 * n;
    /* Stack-local scratch (thread-safe; no static storage). */
    double r [256], p [256], Ap [256];
    for (int i = 0; i < n; i++) {
        if (!ml_isfinite (b [i])) return ML_ERR_NAN_INPUT;
        for (int j = 0; j < n; j++)
            if (!ml_isfinite (ML_TENSOR_AT (A, i, j))) return ML_ERR_NAN_INPUT;
        x [i] = 0.0;
        r [i] = b [i];
        p [i] = b [i];
    }
    long double rs = 0;
    for (int i = 0; i < n; i++) rs += (long double) r [i] * r [i];
    if (!(rs > 0)) return ML_SUCCESS;
    long double bnorm = __builtin_sqrtl (rs);
    /* DESPOT-AUDIT: relative tolerance ||r|| <= tol*||b|| (was absolute
     * ||r|| <= tol, meaningless across scales). tol is now scale-aware. */
    for (int it = 0; it < max_iter; it++) {
        for (int i = 0; i < n; i++) {
            long double s = 0;
            for (int j = 0; j < n; j++) s += (long double) ML_TENSOR_AT (A, i, j) * p [j];
            Ap [i] = (double) s;
        }
        long double pAp = 0;
        for (int i = 0; i < n; i++) pAp += (long double) p [i] * Ap [i];
        if (!(pAp > 0)) return ML_ERR_SINGULAR;
        double alpha = (double) (rs / pAp);
        long double rsn = 0;
        for (int i = 0; i < n; i++) {
            x [i] += alpha * p [i];
            r [i] -= alpha * Ap [i];
            rsn += (long double) r [i] * r [i];
        }
        if (__builtin_sqrtl (rsn) <= (long double) tol * bnorm) return ML_SUCCESS;
        double beta = (double) (rsn / rs);
        for (int i = 0; i < n; i++) p [i] = r [i] + beta * p [i];
        rs = rsn;
    }
    /* R5 despot: exhaustion without meeting tol is CONVERGENCE, not SINGULAR
     * (true breakdown above keeps SINGULAR). */
    return ML_ERR_CONVERGENCE;
}
ML_API ml_status_t ml_gmres_solve (ml_tensor_view_t A, const double *b, double *x, int restart, int max_iter, double tol) {
    int n = A.rows;
    if (!A.data || !b || !x || n <= 0 || n > 64) return ML_ERR_INVALID_ARG;
    if (A.cols != n) return ML_ERR_INVALID_ARG;
    if (restart <= 0) restart = n < 20 ? n : 20;
    if (restart > 32) restart = 32;
    if (!(tol > 0.0) || !ml_isfinite (tol)) return ML_ERR_INVALID_ARG;
    if (max_iter <= 0) max_iter = 4 * n;
    for (int i = 0; i < n; i++) {
        if (!ml_isfinite (b [i])) return ML_ERR_NAN_INPUT;
        for (int j = 0; j < n; j++)
            if (!ml_isfinite (ML_TENSOR_AT (A, i, j))) return ML_ERR_NAN_INPUT;
        x [i] = 0.0;
    }
    long double bnorm0 = 0;
    for (int i = 0; i < n; i++) bnorm0 += (long double) b [i] * b [i];
    bnorm0 = __builtin_sqrtl (bnorm0);
    if (!(bnorm0 > 0)) return ML_SUCCESS;
    /* Stack-local scratch (thread-safe; no static storage). */
    double r [64], V [33][64], H [33][32], cs [32], sn [32], g [33], ys [32];
    for (int outer = 0; outer < max_iter; outer++) {
        for (int i = 0; i < n; i++) {
            long double s = (long double) b [i];
            for (int j = 0; j < n; j++) s -= (long double) ML_TENSOR_AT (A, i, j) * x [j];
            r [i] = (double) s;
        }
        long double beta = 0;
        for (int i = 0; i < n; i++) beta += (long double) r [i] * r [i];
        beta = __builtin_sqrtl (beta);
        if (beta <= (long double) tol * bnorm0) return ML_SUCCESS;
        for (int i = 0; i < n; i++) V [0][i] = (double) ((long double) r [i] / beta);
        for (int i = 0; i <= restart; i++) { g [i] = 0; }
        g [0] = (double) beta;
        int j = 0;
        for (; j < restart; j++) {
            double w [64];
            for (int i = 0; i < n; i++) {
                long double s = 0;
                for (int k = 0; k < n; k++) s += (long double) ML_TENSOR_AT (A, i, k) * V [j][k];
                w [i] = (double) s;
            }
            for (int i = 0; i <= j; i++) {
                long double h = 0;
                for (int k = 0; k < n; k++) h += (long double) w [k] * V [i][k];
                H [i][j] = (double) h;
                for (int k = 0; k < n; k++) w [k] -= H [i][j] * V [i][k];
            }
            long double nw = 0;
            for (int k = 0; k < n; k++) nw += (long double) w [k] * w [k];
            H [j + 1][j] = (double) __builtin_sqrtl (nw);
            if (H [j + 1][j] != 0)
                for (int k = 0; k < n; k++) V [j + 1][k] = (double) ((long double) w [k] / H [j + 1][j]);
            /* Givens */
            for (int i = 0; i < j; i++) {
                double t1 = cs [i] * H [i][j] + sn [i] * H [i + 1][j];
                double t2 = -sn [i] * H [i][j] + cs [i] * H [i + 1][j];
                H [i][j] = t1;
                H [i + 1][j] = t2;
            }
            double den = ml_hypot (H [j][j], H [j + 1][j]);
            if (den == 0) break;
            cs [j] = H [j][j] / den;
            sn [j] = H [j + 1][j] / den;
            H [j][j] = den;
            H [j + 1][j] = 0;
            g [j + 1] = -sn [j] * g [j];
            g [j] *= cs [j];
            if (ml_fabs (g [j + 1]) <= tol * (double) bnorm0) {
                j++;
                break;
            }
        }
        for (int i = j - 1; i >= 0; i--) {
            ys [i] = g [i];
            for (int k = i + 1; k < j; k++) ys [i] -= H [i][k] * ys [k];
            if (H [i][i] == 0) return ML_ERR_SINGULAR;
            ys [i] /= H [i][i];
        }
        for (int i = 0; i < j; i++)
            for (int k = 0; k < n; k++) x [k] += ys [i] * V [i][k];
        if (ml_fabs (g [j]) <= tol * (double) bnorm0) return ML_SUCCESS;
    }
    return ML_ERR_CONVERGENCE;
}
ML_API ml_status_t ml_power_iter (ml_tensor_view_t A, double *lambda, double *vec, int max_iter, double tol) {
    int n = A.rows;
    if (!A.data || !lambda || n <= 0 || n > 128) return ML_ERR_INVALID_ARG;
    if (A.cols != n) return ML_ERR_INVALID_ARG;
    if (max_iter <= 0) max_iter = 1000;
    if (!(tol > 0.0) || !ml_isfinite (tol)) tol = 1e-10;
    /* Stack-local scratch (thread-safe; no static storage). */
    double v [128], w [128];
    for (int i = 0; i < n; i++) v [i] = 1.0 / __builtin_sqrt ((double) n);
    double lam = 0;
    for (int it = 0; it < max_iter; it++) {
        for (int i = 0; i < n; i++) {
            long double s = 0;
            for (int j = 0; j < n; j++) s += (long double) ML_TENSOR_AT (A, i, j) * v [j];
            w [i] = (double) s;
        }
        long double nw = 0;
        for (int i = 0; i < n; i++) nw += (long double) w [i] * w [i];
        nw = __builtin_sqrtl (nw);
        if (!(nw > 0)) return ML_ERR_SINGULAR;
        for (int i = 0; i < n; i++) w [i] = (double) ((long double) w [i] / nw);
        long double rq = 0;
        for (int i = 0; i < n; i++) {
            long double s = 0;
            for (int j = 0; j < n; j++) s += (long double) ML_TENSOR_AT (A, i, j) * w [j];
            rq += (long double) w [i] * s;
        }
        if (ml_fabs ((double) rq - lam) <= tol * (1.0 + ml_fabs ((double) rq))) {
            *lambda = (double) rq;
            if (vec)
                for (int i = 0; i < n; i++) vec [i] = w [i];
            return ML_SUCCESS;
        }
        lam = (double) rq;
        for (int i = 0; i < n; i++) v [i] = w [i];
    }
    return ML_ERR_CONVERGENCE;
}
ML_API ml_status_t ml_svd_jacobi (ml_tensor_view_t A, double *svals, ml_tensor_view_t U, ml_tensor_view_t Vt, int max_sweeps) {
    int m = A.rows, n = A.cols;
    if (!A.data || !svals || m <= 0 || n <= 0 || m > 32 || n > 32) return ML_ERR_INVALID_ARG;
    if (max_sweeps <= 0) max_sweeps = 20;
    /* Stack-local scratch (thread-safe; no static storage). One-sided Jacobi
     * orthogonalizes B columns; V accumulates rotations. UU removed
     * (despot audit): it was dead storage never used for U output. */
    double B [32][32], VV [32][32];
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++) {
            double v = ML_TENSOR_AT (A, i, j);
            if (!ml_isfinite (v)) return ML_ERR_NAN_INPUT;
            B [i][j] = v;
        }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) VV [i][j] = (i == j) ? 1.0 : 0.0;
    int p = m < n ? m : n;
    /* DESPOT-AUDIT: scale-invariant stop. Previous absolute 1e-30 failed
     * for |A|>>1/<<1. Now off <= (eps*||A||_F)^2 with eps=2^-52. */
    long double frob = 0;
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++) frob += (long double) B [i][j] * B [i][j];
    long double off_tol = 2.22044604925031308085e-16L;
    off_tol = off_tol * off_tol * frob;
    for (int sw = 0; sw < max_sweeps; sw++) {
        long double off = 0;
        for (int i = 0; i < p; i++)
            for (int j = i + 1; j < n && j < 32; j++) {
                long double d = 0;
                for (int k = 0; k < m; k++) d += (long double) B [k][i] * B [k][j];
                off += d * d;
            }
        if (!(off > off_tol)) break;
        for (int i = 0; i < n - 1; i++)
            for (int j = i + 1; j < n; j++) {
                long double aii = 0, ajj = 0, aij = 0;
                for (int k = 0; k < m; k++) {
                    aii += (long double) B [k][i] * B [k][i];
                    ajj += (long double) B [k][j] * B [k][j];
                    aij += (long double) B [k][i] * B [k][j];
                }
                if ((double) aij == 0.0) continue;
                double tau = (double) ((ajj - aii) / (2.0L * aij));
                double t = ((tau >= 0) ? 1.0 : -1.0) / (ml_fabs (tau) + __builtin_sqrt (tau * tau + 1.0));
                double c = 1.0 / __builtin_sqrt (t * t + 1.0), s = t * c;
                for (int k = 0; k < m; k++) {
                    double bi = B [k][i], bj = B [k][j];
                    B [k][i] = c * bi - s * bj;
                    B [k][j] = s * bi + c * bj;
                }
                for (int k = 0; k < n; k++) {
                    double bi = VV [k][i], bj = VV [k][j];
                    VV [k][i] = c * bi - s * bj;
                    VV [k][j] = s * bi + c * bj;
                }
            }
    }
    for (int j = 0; j < n; j++) {
        long double s = 0;
        for (int k = 0; k < m; k++) s += (long double) B [k][j] * B [k][j];
        svals [j] = (double) __builtin_sqrtl (s);
        if (svals [j] != 0 && U.data)
            for (int k = 0; k < m; k++) ML_TENSOR_AT (U, k, j) = (double) ((long double) B [k][j] / svals [j]);
    }
    if (Vt.data) {
        if (Vt.rows != n || Vt.cols != n) return ML_ERR_INVALID_ARG;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) ML_TENSOR_AT (Vt, i, j) = VV [j][i];
    }
    /* Sort singular values descending, permuting U columns and Vt rows
     * consistently via the same pairwise swaps (permutation tracking).
     * Previously only svals were swapped, leaving U/Vt mismatched. */
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (svals [j] > svals [i]) {
                double t = svals [i];
                svals [i] = svals [j];
                svals [j] = t;
                if (U.data) {
                    for (int k = 0; k < m; k++) {
                        double tu = ML_TENSOR_AT (U, k, i);
                        ML_TENSOR_AT (U, k, i) = ML_TENSOR_AT (U, k, j);
                        ML_TENSOR_AT (U, k, j) = tu;
                    }
                }
                if (Vt.data) {
                    for (int k = 0; k < n; k++) {
                        double tv = ML_TENSOR_AT (Vt, i, k);
                        ML_TENSOR_AT (Vt, i, k) = ML_TENSOR_AT (Vt, j, k);
                        ML_TENSOR_AT (Vt, j, k) = tv;
                    }
                }
            }
        }
    }
    return ML_SUCCESS;
}
ML_API ml_status_t ml_qr_iter_eig (ml_tensor_view_t A, double *evals_re, double *evals_im, int max_iter, double tol) {
    int n = A.rows;
    if (!A.data || !evals_re || !evals_im || n <= 0 || n > 16) return ML_ERR_INVALID_ARG;
    if (A.cols != n) return ML_ERR_INVALID_ARG;
    if (!(tol > 0.0) || !ml_isfinite (tol)) tol = 1e-10;
    if (max_iter <= 0) max_iter = 500;
    /* Stack-local scratch (thread-safe; no static storage). */
    double H [16][16];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            double v = ML_TENSOR_AT (A, i, j);
            if (!ml_isfinite (v)) return ML_ERR_NAN_INPUT;
            H [i][j] = v;
        }
    /* R6 despot: Householder reduction to upper Hessenberg form first.
     * The old code swept Givens rotations on the FULL matrix with each
     * right-rotation applied immediately, nesting similarities instead of
     * forming RQ: it converged on tridiagonal/diagonal inputs but stalled
     * (~20% failure) on general SPD, including 2x2 rotation blocks already
     * in real Schur form. Francis QR below is the textbook algorithm. */
    for (int m = 0; m < n - 1; m++) {
        long double x0 = (long double) H [m + 1][m];
        long double sigma = 0.0L;
        long double mx = fabsl (x0);
        int k;
        for (k = m + 2; k < n; k++) {
            long double v = fabsl ((long double) H [k][m]);
            if (v > mx) mx = v;
        }
        if (mx == 0.0L) continue;
        for (k = m + 1; k < n; k++) {
            long double v = (long double) H [k][m] / mx;
            sigma += v * v;
        }
        sigma = mx * __builtin_sqrtl (sigma);
        if (sigma == 0.0L) continue;
        {
            /* Reflector P = I - beta u u', u[0] = 1, alpha = -sgn(x0)*|x|. */
            long double u [16];
            long double alpha = (x0 <= 0.0L) ? sigma : -sigma;
            long double v0 = x0 - alpha;
            long double beta;
            int len = n - m - 1;
            if (v0 == 0.0L) continue;
            for (k = 0; k < len; k++) u [k] = (long double) H [m + 1 + k][m] / v0;
            u [0] = 1.0L;
            /* R6 despot: beta is -v0/alpha for the v0-normalized u
             * (2/||u||^2 with ||u||^2 = -2 alpha/v0); the old -1/(alpha v0)
             * was off by v0^2 and silently broke similarity (trace drift). */
            beta = -v0 / alpha;
            for (int j = m; j < n; j++) {
                long double dot = 0.0L;
                for (k = 0; k < len; k++) dot += u [k] * (long double) H [m + 1 + k][j];
                dot *= beta;
                for (k = 0; k < len; k++) H [m + 1 + k][j] = (double) ((long double) H [m + 1 + k][j] - dot * u [k]);
            }
            for (int i2 = 0; i2 < n; i2++) {
                long double dot = 0.0L;
                for (k = 0; k < len; k++) dot += (long double) H [i2][m + 1 + k] * u [k];
                dot *= beta;
                for (k = 0; k < len; k++) H [i2][m + 1 + k] = (double) ((long double) H [i2][m + 1 + k] - dot * u [k]);
            }
        }
    }
    {
        double out_re [16];
        double out_im [16];
        int nout = 0;
        int hi = n - 1;
        int sweeps = 0;
        int stall = 0;
        while (hi >= 0) {
            int lo;
            if (hi == 0) {
                out_re [nout] = H [0][0];
                out_im [nout] = 0.0;
                nout++;
                break;
            }
            if (!(ml_fabs (H [hi][hi - 1]) > tol)) {
                out_re [nout] = H [hi][hi];
                out_im [nout] = 0.0;
                nout++;
                hi--;
                continue;
            }
            lo = hi - 1;
            while (lo > 0 && ml_fabs (H [lo][lo - 1]) > tol) lo--;
            if (hi == lo + 1) {
                /* Decoupled 2x2: solve directly (real or complex pair). */
                long double a = (long double) H [lo][lo];
                long double b = (long double) H [lo][hi];
                long double c = (long double) H [hi][lo];
                long double d = (long double) H [hi][hi];
                long double tr = a + d;
                long double det = a * d - b * c;
                long double disc = tr * tr - 4.0L * det;
                if (disc < 0.0L) {
                    double s = (double) (__builtin_sqrtl (-disc) / 2.0L);
                    out_re [nout] = (double) (tr / 2.0L);
                    out_im [nout] = s;
                    nout++;
                    out_re [nout] = (double) (tr / 2.0L);
                    out_im [nout] = -s;
                    nout++;
                } else {
                    double s = (double) (__builtin_sqrtl (disc) / 2.0L);
                    out_re [nout] = (double) (tr / 2.0L) + s;
                    out_im [nout] = 0.0;
                    nout++;
                    out_re [nout] = (double) (tr / 2.0L) - s;
                    out_im [nout] = 0.0;
                    nout++;
                }
                hi = lo - 1;
                continue;
            }
            /* Explicit shifted QR on [lo..hi] (R6 despot).
             * The Francis bulge chase it replaces left (k+3,k+1)-type fill
             * below the band that no chase step ever cleaned (~50% stall
             * on general 4x4+). Here Q = orthog(M) via MGS and H <- Q'HQ is
             * a full-block similarity by construction: no band to maintain,
             * nothing to chase. Normal path is the double shift
             * M = H^2-tr H+det I (bottom 2x2; handles complex pairs that
             * stalled single-shift Wilkinson). Exceptional path (stall>10)
             * is a SINGLE ad-hoc shift M = H-mu I: it breaks the exact
             * period-2 cycles where |p(nd)| == 1 on the whole spectrum
             * renders every double shift spectrally invisible. */
            {
                double tr = 0.0;
                double det = 0.0;
                double mu = 0.0;
                int single = (stall > 10);
                double M [16][16];
                double Q [16][16];
                double T [16][16];
                int m = hi - lo + 1;
                int i;
                int j;
                int k;
                if (single) {
                    stall = 0;
                    /* Sweep blocks reaching here are >= 3x3 (2x2 solves
                     * directly), so H[hi-1][hi-2] is in-block. */
                    mu = ml_fabs (H [hi][hi - 1]) + ml_fabs (H [hi - 1][hi - 2]);
                    if (!(mu > 0.0) || !ml_isfinite (mu)) mu = ml_fabs (H [hi][hi]) + 1.0;
                } else {
                    long double a = (long double) H [hi - 1][hi - 1];
                    long double b = (long double) H [hi - 1][hi];
                    long double c = (long double) H [hi][hi - 1];
                    long double d = (long double) H [hi][hi];
                    tr = (double) (a + d);
                    det = (double) (a * d - b * c);
                }
                if (!ml_isfinite (tr) || !ml_isfinite (det) || !ml_isfinite (mu)) {
                    for (k = 0; k < n; k++) {
                        evals_re [k] = ml_make_nan ();
                        evals_im [k] = ml_make_nan ();
                    }
                    return ML_ERR_SINGULAR;
                }
                /* M = shift polynomial on the block, LD accumulation. */
                for (i = 0; i < m; i++) {
                    for (j = 0; j < m; j++) {
                        long double s;
                        if (single) {
                            s = (long double) H [lo + i][lo + j];
                            if (i == j) s -= (long double) mu;
                        } else {
                            s = 0.0L;
                            for (k = 0; k < m; k++) s += (long double) H [lo + i][lo + k] * (long double) H [lo + k][lo + j];
                            s -= (long double) tr * (long double) H [lo + i][lo + j];
                            if (i == j) s += (long double) det;
                        }
                        M [i][j] = (double) s;
                        if (!ml_isfinite (M [i][j])) {
                            for (k = 0; k < n; k++) {
                                evals_re [k] = ml_make_nan ();
                                evals_im [k] = ml_make_nan ();
                            }
                            return ML_ERR_OVERFLOW;
                        }
                    }
                }
                /* MGS orthogonalization of M's columns, done TWICE (R6).
                 * Classical MGS loses orthogonality (~eps*cond(M)) precisely
                 * when p(H) goes singular near convergence; the resulting
                 * non-orthogonal "similarities" then inject drift (~1e-8)
                 * and stall deflation for hundreds of sweeps. A second
                 * pass (Giraud-Langou) restores orthogonality to ~eps.
                 * Dependent columns get a re-orthogonalized unit fallback,
                 * never divide-by-zero. R factor unneeded. */
                for (i = 0; i < m; i++) {
                    for (j = 0; j < m; j++) {
                        Q [i][j] = 0.0;
                    }
                }
                for (j = 0; j < m; j++) {
                    long double v [16];
                    int pass;
                    for (i = 0; i < m; i++) v [i] = (long double) M [i][j];
                    for (pass = 0; pass < 2; pass++) {
                        for (k = 0; k < j; k++) {
                            long double dot = 0.0L;
                            for (i = 0; i < m; i++) dot += (long double) Q [i][k] * v [i];
                            for (i = 0; i < m; i++) v [i] -= dot * (long double) Q [i][k];
                        }
                    }
                    {
                        long double nn = 0.0L;
                        for (i = 0; i < m; i++) nn += v [i] * v [i];
                        nn = __builtin_sqrtl (nn);
                        if (!(nn > 1e-300L) || !ml_isfinite ((double) nn)) {
                            int t;
                            int ok = 0;
                            for (t = 0; t < m && !ok; t++) {
                                long double w [16] = {0.0L};
                                long double wn;
                                int q;
                                w [t] = 1.0L;
                                for (q = 0; q < j; q++) {
                                    long double dot = 0.0L;
                                    for (i = 0; i < m; i++) dot += (long double) Q [i][q] * w [i];
                                    for (i = 0; i < m; i++) w [i] -= dot * (long double) Q [i][q];
                                }
                                wn = 0.0L;
                                for (i = 0; i < m; i++) wn += w [i] * w [i];
                                wn = __builtin_sqrtl (wn);
                                if (wn > 0.5L && ml_isfinite ((double) wn)) {
                                    for (i = 0; i < m; i++) v [i] = w [i] / wn;
                                    nn = 1.0L;
                                    ok = 1;
                                }
                            }
                            if (!ok) {
                                for (k = 0; k < n; k++) {
                                    evals_re [k] = ml_make_nan ();
                                    evals_im [k] = ml_make_nan ();
                                }
                                return ML_ERR_SINGULAR;
                            }
                        }
                        for (i = 0; i < m; i++) Q [i][j] = (double) (v [i] / nn);
                    }
                }
                /* H <- Q' H Q on the block. */
                for (i = 0; i < m; i++) {
                    for (j = 0; j < m; j++) {
                        long double s = 0.0L;
                        for (k = 0; k < m; k++) s += (long double) Q [k][i] * (long double) H [lo + k][lo + j];
                        T [i][j] = (double) s;
                    }
                }
                for (i = 0; i < m; i++) {
                    for (j = 0; j < m; j++) {
                        long double s = 0.0L;
                        for (k = 0; k < m; k++) s += T [i][k] * (long double) Q [k][j];
                        if (!ml_isfinite ((double) s)) {
                            for (k = 0; k < n; k++) {
                                evals_re [k] = ml_make_nan ();
                                evals_im [k] = ml_make_nan ();
                            }
                            return ML_ERR_OVERFLOW;
                        }
                        H [lo + i][lo + j] = (double) s;
                    }
                }
            }
            sweeps++;
            stall++;
            if (sweeps > max_iter) {
                for (int k = 0; k < n; k++) {
                    evals_re [k] = ml_make_nan ();
                    evals_im [k] = ml_make_nan ();
                }
                return ML_ERR_CONVERGENCE;
            }
        }
        for (int k = 0; k < n; k++) {
            evals_re [k] = (k < nout) ? out_re [k] : ml_make_nan ();
            evals_im [k] = (k < nout) ? out_im [k] : ml_make_nan ();
        }
        return ML_SUCCESS;
    }
}
