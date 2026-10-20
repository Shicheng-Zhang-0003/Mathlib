#include "ml_compiler.h"
#include "ml_linalg_extra.h"
#include "ml_exp_log.h"
#include <stddef.h>
/* v12R2 linear-algebra depth pack: Hessenberg, SVD, pseudoinverse,
 * condition number, rank, iterative solvers and matrix functions.
 * Stateless; all scratch comes from the caller's ml_workspace_t. */

static void mm_identity (double *M, int n) {
    int i, j;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            M [(size_t) i * (size_t) n + (size_t) j] = (i == j) ? 1.0 : 0.0;
}
/* C = A * B  (row-major n x n). C must not alias A or B. */
static void mm_mul (double *C, const double *A, const double *B, int n) {
    int i, j, k;
    size_t sn = (size_t) n;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            double s = 0.0;
            for (k = 0; k < n; k++)
                s += A [(size_t) i * sn + (size_t) k] * B [(size_t) k * sn + (size_t) j];
            C [(size_t) i * sn + (size_t) j] = s;
        }
    }
}
static double mm_norm1 (const double *A, int n) {
    double mx = 0.0;
    int i, j;
    size_t sn = (size_t) n;
    for (j = 0; j < n; j++) {
        double s = 0.0;
        for (i = 0; i < n; i++) s += ml_fabs (A [(size_t) i * sn + (size_t) j]);
        if (s > mx) mx = s;
    }
    return mx;
}
/* Gauss-Jordan inverse with partial pivoting. work is n x n scratch.
 * Returns 0 on success, -1 if numerically singular. */
static int mm_inv (const double *A, double *Inv, int n, double *work) {
    size_t sn = (size_t) n;
    int i, j, col;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            work [(size_t) i * sn + (size_t) j] = A [(size_t) i * sn + (size_t) j];
    mm_identity (Inv, n);
    for (col = 0; col < n; col++) {
        int piv = col;
        double best = ml_fabs (work [(size_t) col * sn + (size_t) col]);
        for (i = col + 1; i < n; i++) {
            double v = ml_fabs (work [(size_t) i * sn + (size_t) col]);
            if (v > best) {
                best = v;
                piv = i;
            }
        }
        if (best == 0.0) return -1;
        if (piv != col) {
            for (j = 0; j < n; j++) {
                double t = work [(size_t) col * sn + (size_t) j];
                work [(size_t) col * sn + (size_t) j] = work [(size_t) piv * sn + (size_t) j];
                work [(size_t) piv * sn + (size_t) j] = t;
                t = Inv [(size_t) col * sn + (size_t) j];
                Inv [(size_t) col * sn + (size_t) j] = Inv [(size_t) piv * sn + (size_t) j];
                Inv [(size_t) piv * sn + (size_t) j] = t;
            }
        }
        {
            double d = work [(size_t) col * sn + (size_t) col];
            for (j = 0; j < n; j++) {
                work [(size_t) col * sn + (size_t) j] /= d;
                Inv [(size_t) col * sn + (size_t) j] /= d;
            }
        }
        for (i = 0; i < n; i++) {
            double f;
            if (i == col) continue;
            f = work [(size_t) i * sn + (size_t) col];
            if (f == 0.0) continue;
            for (j = 0; j < n; j++) {
                work [(size_t) i * sn + (size_t) j] -= f * work [(size_t) col * sn + (size_t) j];
                Inv [(size_t) i * sn + (size_t) j] -= f * Inv [(size_t) col * sn + (size_t) j];
            }
        }
    }
    return 0;
}
/* ------------------------------------------------------------------ SVD --- */

/* One-sided Jacobi SVD on a column-major R x C matrix W (R >= C). On exit W
 * holds U (columns normalized), V (column-major C x C) holds the right
 * singular vectors and S[C] the singular values, sorted descending. */
static void svd_jacobi (double *W, int R, int C, double *V, double *S) {
    int i, p, q, sweep;
    size_t Rz = (size_t) R, Cz = (size_t) C;
    for (i = 0; i < C; i++) {
        int j;
        for (j = 0; j < C; j++)
            V [(size_t) i + Cz * (size_t) j] = (i == j) ? 1.0 : 0.0;
    }
    for (sweep = 0; sweep < 60; sweep++) {
        int changed = 0;
        for (p = 0; p < C - 1; p++) {
            for (q = p + 1; q < C; q++) {
                double alpha = 0.0, beta = 0.0, gamma = 0.0;
                double zeta, t, c, s;
                for (i = 0; i < R; i++) {
                    double ap = W [(size_t) i + Rz * (size_t) p];
                    double aq = W [(size_t) i + Rz * (size_t) q];
                    alpha += ap * ap;
                    beta += aq * aq;
                    gamma += ap * aq;
                }
                if (gamma == 0.0) continue;
                if (ml_fabs (gamma) <= 1e-17 * ml_sqrt (alpha * beta)) continue;
                zeta = (beta - alpha) / (2.0 * gamma);
                t = ((zeta >= 0.0) ? 1.0 : -1.0) / (ml_fabs (zeta) + ml_sqrt (1.0 + zeta * zeta));
                c = 1.0 / ml_sqrt (1.0 + t * t);
                s = c * t;
                for (i = 0; i < R; i++) {
                    double ap = W [(size_t) i + Rz * (size_t) p];
                    double aq = W [(size_t) i + Rz * (size_t) q];
                    W [(size_t) i + Rz * (size_t) p] = c * ap - s * aq;
                    W [(size_t) i + Rz * (size_t) q] = s * ap + c * aq;
                }
                for (i = 0; i < C; i++) {
                    double vp = V [(size_t) i + Cz * (size_t) p];
                    double vq = V [(size_t) i + Cz * (size_t) q];
                    V [(size_t) i + Cz * (size_t) p] = c * vp - s * vq;
                    V [(size_t) i + Cz * (size_t) q] = s * vp + c * vq;
                }
                changed = 1;
            }
        }
        if (!changed) break;
    }
    for (q = 0; q < C; q++) {
        double nrm = 0.0;
        for (i = 0; i < R; i++) {
            double a = W [(size_t) i + Rz * (size_t) q];
            nrm += a * a;
        }
        nrm = ml_sqrt (nrm);
        S [q] = nrm;
        if (nrm > 0.0)
            for (i = 0; i < R; i++) W [(size_t) i + Rz * (size_t) q] /= nrm;
    }
    for (p = 0; p < C - 1; p++) {
        int mx = p;
        for (q = p + 1; q < C; q++)
            if (S [q] > S [mx]) mx = q;
        if (mx != p) {
            double tmp = S [p];
            S [p] = S [mx];
            S [mx] = tmp;
            for (i = 0; i < R; i++) {
                double t = W [(size_t) i + Rz * (size_t) p];
                W [(size_t) i + Rz * (size_t) p] = W [(size_t) i + Rz * (size_t) mx];
                W [(size_t) i + Rz * (size_t) mx] = t;
            }
            for (i = 0; i < C; i++) {
                double t = V [(size_t) i + Cz * (size_t) p];
                V [(size_t) i + Cz * (size_t) p] = V [(size_t) i + Cz * (size_t) mx];
                V [(size_t) i + Cz * (size_t) mx] = t;
            }
        }
    }
}
ML_API ml_status_t ml_svd (ml_tensor_view_t A, ml_tensor_view_t U, double *S,
                           ml_tensor_view_t Vt, ml_workspace_t *ws) {
    int m = A.rows, n = A.cols, k, R, i, j;
    size_t need;
    double *W, *V;
    if (ML_UNLIKELY (m <= 0 || n <= 0 || !A.data || !S || !ws)) return ML_ERR_INVALID_ARG;
    k = (m < n) ? m : n;
    R = (m > n) ? m : n;
    if (U.data && (U.rows != m || U.cols != k)) return ML_ERR_INVALID_ARG;
    if (Vt.data && (Vt.rows != k || Vt.cols != n)) return ML_ERR_INVALID_ARG;
    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            if (ML_UNLIKELY (!ml_isfinite (ML_TENSOR_AT (A, i, j)))) return ML_ERR_NAN_INPUT;
    need = ((size_t) R * (size_t) k + (size_t) k * (size_t) k) * sizeof (double);
    if (ML_UNLIKELY (ws->capacity < need)) return ML_ERR_WORKSPACE;
    W = (double *) ml_workspace_alloc (ws, (size_t) R * (size_t) k * sizeof (double));
    V = (double *) ml_workspace_alloc (ws, (size_t) k * (size_t) k * sizeof (double));
    if (ML_UNLIKELY (!W || !V)) return ML_ERR_WORKSPACE;
    if (m >= n) {
        for (j = 0; j < n; j++)
            for (i = 0; i < m; i++)
                W [(size_t) i + (size_t) m * (size_t) j] = ML_TENSOR_AT (A, i, j);
    } else {
        for (j = 0; j < m; j++)
            for (i = 0; i < n; i++)
                W [(size_t) i + (size_t) n * (size_t) j] = ML_TENSOR_AT (A, j, i);
    }
    svd_jacobi (W, R, k, V, S);
    if (m >= n) {
        if (U.data)
            for (j = 0; j < k; j++)
                for (i = 0; i < m; i++)
                    ML_TENSOR_AT (U, i, j) = W [(size_t) i + (size_t) m * (size_t) j];
        if (Vt.data)
            for (j = 0; j < k; j++)
                for (i = 0; i < n; i++)
                    ML_TENSOR_AT (Vt, j, i) = V [(size_t) i + (size_t) n * (size_t) j];
    } else {
        if (U.data)
            for (j = 0; j < k; j++)
                for (i = 0; i < m; i++)
                    ML_TENSOR_AT (U, i, j) = V [(size_t) i + (size_t) m * (size_t) j];
        if (Vt.data)
            for (j = 0; j < n; j++)
                for (i = 0; i < k; i++)
                    ML_TENSOR_AT (Vt, i, j) = W [(size_t) j + (size_t) n * (size_t) i];
    }
    return ML_SUCCESS;
}
ML_API ml_status_t ml_cond2 (ml_tensor_view_t A, double *out_cond,
                             ml_workspace_t *ws) {
    int m = A.rows, n = A.cols, k;
    double *S;
    ml_status_t st;
    if (ML_UNLIKELY (!out_cond)) return ML_ERR_INVALID_ARG;
    k = (m < n) ? m : n;
    if (ML_UNLIKELY (k <= 0 || !ws)) return ML_ERR_INVALID_ARG;
    S = (double *) ml_workspace_alloc (ws, (size_t) k * sizeof (double));
    if (ML_UNLIKELY (!S)) return ML_ERR_WORKSPACE;
    st = ml_svd (A, ml_tensor_view (NULL, 0, 0), S, ml_tensor_view (NULL, 0, 0), ws);
    if (st != ML_SUCCESS) return st;
    if (S [k - 1] <= 0.0) {
        *out_cond = ml_make_inf (0);
        return ML_SUCCESS;
    }
    *out_cond = S [0] / S [k - 1];
    return ML_SUCCESS;
}
ML_API ml_status_t ml_matrix_rank (ml_tensor_view_t A, double tol, int *rank,
                                   ml_workspace_t *ws) {
    int m = A.rows, n = A.cols, k, i;
    double *S;
    ml_status_t st;
    if (ML_UNLIKELY (!rank)) return ML_ERR_INVALID_ARG;
    k = (m < n) ? m : n;
    if (ML_UNLIKELY (k <= 0 || !ws)) return ML_ERR_INVALID_ARG;
    S = (double *) ml_workspace_alloc (ws, (size_t) k * sizeof (double));
    if (ML_UNLIKELY (!S)) return ML_ERR_WORKSPACE;
    st = ml_svd (A, ml_tensor_view (NULL, 0, 0), S, ml_tensor_view (NULL, 0, 0), ws);
    if (st != ML_SUCCESS) return st;
    *rank = 0;
    for (i = 0; i < k; i++)
        if (S [i] > tol) (*rank)++;
    return ML_SUCCESS;
}
ML_API ml_status_t ml_pinv (ml_tensor_view_t A, ml_tensor_view_t Apinv,
                            double tol, ml_workspace_t *ws) {
    int m = A.rows, n = A.cols, k, R, i, j, l;
    size_t need;
    double *W, *V, *S;
    if (ML_UNLIKELY (m <= 0 || n <= 0 || !A.data || !Apinv.data || !ws)) return ML_ERR_INVALID_ARG;
    if (ML_UNLIKELY (Apinv.rows != n || Apinv.cols != m)) return ML_ERR_INVALID_ARG;
    k = (m < n) ? m : n;
    R = (m > n) ? m : n;
    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            if (ML_UNLIKELY (!ml_isfinite (ML_TENSOR_AT (A, i, j)))) return ML_ERR_NAN_INPUT;
    need = ((size_t) R * (size_t) k + (size_t) k * (size_t) k + (size_t) k) * sizeof (double);
    if (ML_UNLIKELY (ws->capacity < need)) return ML_ERR_WORKSPACE;
    W = (double *) ml_workspace_alloc (ws, (size_t) R * (size_t) k * sizeof (double));
    V = (double *) ml_workspace_alloc (ws, (size_t) k * (size_t) k * sizeof (double));
    S = (double *) ml_workspace_alloc (ws, (size_t) k * sizeof (double));
    if (ML_UNLIKELY (!W || !V || !S)) return ML_ERR_WORKSPACE;
    if (m >= n) {
        for (j = 0; j < n; j++)
            for (i = 0; i < m; i++)
                W [(size_t) i + (size_t) m * (size_t) j] = ML_TENSOR_AT (A, i, j);
    } else {
        for (j = 0; j < m; j++)
            for (i = 0; i < n; i++)
                W [(size_t) i + (size_t) n * (size_t) j] = ML_TENSOR_AT (A, j, i);
    }
    svd_jacobi (W, R, k, V, S);
    /* pinv(r,c) = sum_l V(r,l) * (1/S_l) * U(c,l), where V(l,r) uses Vt. */
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            double s = 0.0;
            for (l = 0; l < k; l++) {
                double sv = S [l];
                double vl = (m >= n) ? V [(size_t) i + (size_t) n * (size_t) l]
                                     : W [(size_t) i + (size_t) n * (size_t) l];
                double uc = (m >= n) ? W [(size_t) j + (size_t) m * (size_t) l]
                                     : V [(size_t) j + (size_t) m * (size_t) l];
                if (sv > tol) s += vl * uc / sv;
            }
            ML_TENSOR_AT (Apinv, i, j) = s;
        }
    }
    return ML_SUCCESS;
}
/* ---------------------------------------------------------- Hessenberg --- */

ML_API ml_status_t ml_hessenberg (ml_tensor_view_t A, ml_tensor_view_t H,
                                  ml_tensor_view_t Q) {
    int n = A.rows, k, i, j;
    double *v, vstack [256];
    if (ML_UNLIKELY (n <= 0 || !A.data || !H.data)) return ML_ERR_INVALID_ARG;
    if (ML_UNLIKELY (A.cols != n || H.rows != n || H.cols != n)) return ML_ERR_INVALID_ARG;
    if (Q.data && (Q.rows != n || Q.cols != n)) return ML_ERR_INVALID_ARG;
    if (ML_UNLIKELY (n > 256)) return ML_ERR_WORKSPACE;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++) {
            double a = ML_TENSOR_AT (A, i, j);
            if (ML_UNLIKELY (!ml_isfinite (a))) return ML_ERR_NAN_INPUT;
            ML_TENSOR_AT (H, i, j) = a;
        }
    if (Q.data) mm_identity (Q.data, n);
    v = vstack;
    for (k = 0; k < n - 2; k++) {
        int L = n - k - 1;
        double nx = 0.0, alpha, nv2 = 0.0;
        for (i = 0; i < n; i++) v [i] = 0.0;
        for (i = 0; i < L; i++) {
            double x = ML_TENSOR_AT (H, k + 1 + i, k);
            v [k + 1 + i] = x;
            nx += x * x;
        }
        nx = ml_sqrt (nx);
        if (nx == 0.0) continue;
        alpha = ((ML_TENSOR_AT (H, k + 1, k) >= 0.0) ? -nx : nx);
        v [k + 1] -= alpha;
        for (i = 0; i < L; i++) nv2 += v [k + 1 + i] * v [k + 1 + i];
        if (nv2 == 0.0) continue;
        {
            double inv2 = 2.0 / nv2;
            for (j = k; j < n; j++) {
                double w = 0.0;
                for (i = 0; i < L; i++) w += v [k + 1 + i] * ML_TENSOR_AT (H, k + 1 + i, j);
                w *= inv2;
                for (i = 0; i < L; i++)
                    ML_TENSOR_AT (H, k + 1 + i, j) -= v [k + 1 + i] * w;
            }
            for (i = 0; i < n; i++) {
                double u = 0.0;
                for (j = 0; j < L; j++) u += ML_TENSOR_AT (H, i, k + 1 + j) * v [k + 1 + j];
                u *= inv2;
                for (j = 0; j < L; j++)
                    ML_TENSOR_AT (H, i, k + 1 + j) -= u * v [k + 1 + j];
            }
            if (Q.data) {
                for (i = 0; i < n; i++) {
                    double qv = 0.0;
                    for (j = 0; j < L; j++) qv += ML_TENSOR_AT (Q, i, k + 1 + j) * v [k + 1 + j];
                    qv *= inv2;
                    for (j = 0; j < L; j++)
                        ML_TENSOR_AT (Q, i, k + 1 + j) -= qv * v [k + 1 + j];
                }
            }
        }
        for (i = k + 2; i < n; i++) ML_TENSOR_AT (H, i, k) = 0.0;
        ML_TENSOR_AT (H, k + 1, k) = alpha;
    }
    return ML_SUCCESS;
}
/* ---------------------------------------------------- Iterative solvers --- */

ML_API ml_status_t ml_gauss_seidel (ml_tensor_view_t A, const double *b,
                                    double *x, int max_iter, double tol,
                                    int *iters) {
    int n = A.rows, i, j, it;
    if (ML_UNLIKELY (n <= 0 || !A.data || !b || !x)) return ML_ERR_INVALID_ARG;
    if (ML_UNLIKELY (A.cols != n)) return ML_ERR_INVALID_ARG;
    if (max_iter <= 0) max_iter = 1000;
    if (tol <= 0.0) tol = 1e-12;
    for (it = 0; it < max_iter; it++) {
        double diff = 0.0;
        for (i = 0; i < n; i++) {
            double s = b [i], xold = x [i], diag = ML_TENSOR_AT (A, i, i);
            if (diag == 0.0) return ML_ERR_SINGULAR;
            for (j = 0; j < n; j++)
                if (j != i) s -= ML_TENSOR_AT (A, i, j) * x [j];
            x [i] = s / diag;
            {
                double d = ml_fabs (x [i] - xold);
                if (d > diff) diff = d;
            }
        }
        if (diff < tol) {
            if (iters) *iters = it + 1;
            return ML_SUCCESS;
        }
    }
    if (iters) *iters = max_iter;
    return ML_ERR_CONVERGENCE;
}
ML_API ml_status_t ml_sor (ml_tensor_view_t A, const double *b, double *x,
                           double omega, int max_iter, double tol, int *iters) {
    int n = A.rows, i, j, it;
    if (ML_UNLIKELY (n <= 0 || !A.data || !b || !x)) return ML_ERR_INVALID_ARG;
    if (ML_UNLIKELY (A.cols != n)) return ML_ERR_INVALID_ARG;
    if (omega <= 0.0 || omega >= 2.0) return ML_ERR_INVALID_ARG;
    if (max_iter <= 0) max_iter = 1000;
    if (tol <= 0.0) tol = 1e-12;
    for (it = 0; it < max_iter; it++) {
        double diff = 0.0;
        for (i = 0; i < n; i++) {
            double s = b [i], xold = x [i], diag = ML_TENSOR_AT (A, i, i);
            if (diag == 0.0) return ML_ERR_SINGULAR;
            for (j = 0; j < n; j++)
                if (j != i) s -= ML_TENSOR_AT (A, i, j) * x [j];
            x [i] = (1.0 - omega) * xold + omega * (s / diag);
            {
                double d = ml_fabs (x [i] - xold);
                if (d > diff) diff = d;
            }
        }
        if (diff < tol) {
            if (iters) *iters = it + 1;
            return ML_SUCCESS;
        }
    }
    if (iters) *iters = max_iter;
    return ML_ERR_CONVERGENCE;
}
ML_API ml_status_t ml_conjugate_gradient (ml_tensor_view_t A, const double *b,
                                          double *x, int max_iter, double tol,
                                          int *iters) {
    int n = A.rows, i, j, it;
    double rsold;
    double rstack [512], pstack [512], apstack [512];
    double *r, *p, *ap;
    if (ML_UNLIKELY (n <= 0 || !A.data || !b || !x)) return ML_ERR_INVALID_ARG;
    if (ML_UNLIKELY (A.cols != n)) return ML_ERR_INVALID_ARG;
    if (ML_UNLIKELY (n > 512)) return ML_ERR_WORKSPACE;
    if (max_iter <= 0) max_iter = 10 * n;
    if (tol <= 0.0) tol = 1e-12;
    r = rstack;
    p = pstack;
    ap = apstack;
    for (i = 0; i < n; i++) {
        double s = 0.0;
        for (j = 0; j < n; j++) s += ML_TENSOR_AT (A, i, j) * x [j];
        r [i] = b [i] - s;
        p [i] = r [i];
    }
    rsold = 0.0;
    for (i = 0; i < n; i++) rsold += r [i] * r [i];
    if (ml_sqrt (rsold) < tol) {
        if (iters) *iters = 0;
        return ML_SUCCESS;
    }
    for (it = 0; it < max_iter; it++) {
        double pap = 0.0, alpha, rsnew = 0.0;
        for (i = 0; i < n; i++) {
            double s = 0.0;
            for (j = 0; j < n; j++) s += ML_TENSOR_AT (A, i, j) * p [j];
            ap [i] = s;
            pap += p [i] * s;
        }
        if (pap == 0.0) break;
        alpha = rsold / pap;
        for (i = 0; i < n; i++) {
            x [i] += alpha * p [i];
            r [i] -= alpha * ap [i];
            rsnew += r [i] * r [i];
        }
        if (ml_sqrt (rsnew) < tol) {
            if (iters) *iters = it + 1;
            return ML_SUCCESS;
        }
        for (i = 0; i < n; i++) p [i] = r [i] + (rsnew / rsold) * p [i];
        rsold = rsnew;
    }
    if (iters) *iters = max_iter;
    return ML_ERR_CONVERGENCE;
}
/* ------------------------------------------------------ Matrix functions -- */

/* Compute sqrt of A into R using Denman-Beavers. buf holds 6 n*n blocks:
 * Y, Z, Yi, Zi, W, T. A is not modified. */
static int mat_sqrt_db (const double *A, double *R, int n, double *buf) {
    size_t nn = (size_t) n * (size_t) n;
    double *Y = buf, *Z = buf + nn, *Yi = buf + 2 * nn, *Zi = buf + 3 * nn;
    double *W = buf + 4 * nn, *T = buf + 5 * nn;
    int i, it;
    for (i = 0; i < (int) nn; i++) {
        Y [i] = A [i];
        Z [i] = 0.0;
    }
    for (i = 0; i < n; i++) Z [(size_t) i * (size_t) n + (size_t) i] = 1.0;
    for (it = 0; it < 40; it++) {
        double diff = 0.0;
        if (mm_inv (Y, Yi, n, W) != 0) return -1;
        if (mm_inv (Z, Zi, n, W) != 0) return -1;
        for (i = 0; i < (int) nn; i++) {
            T [i] = 0.5 * (Y [i] + Zi [i]);
            W [i] = 0.5 * (Z [i] + Yi [i]);
        }
        for (i = 0; i < (int) nn; i++) {
            double d = ml_fabs (T [i] - Y [i]);
            if (d > diff) diff = d;
        }
        for (i = 0; i < (int) nn; i++) {
            Y [i] = T [i];
            Z [i] = W [i];
        }
        if (diff < 1e-15) break;
    }
    for (i = 0; i < (int) nn; i++) R [i] = Y [i];
    return 0;
}
ML_API ml_status_t ml_matrix_sqrt (ml_tensor_view_t A, ml_tensor_view_t R,
                                   ml_workspace_t *ws) {
    int n = A.rows, i, j;
    size_t nn, need;
    double *buf;
    if (ML_UNLIKELY (n <= 0 || !A.data || !R.data || !ws)) return ML_ERR_INVALID_ARG;
    if (ML_UNLIKELY (A.cols != n || R.rows != n || R.cols != n)) return ML_ERR_INVALID_ARG;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            if (ML_UNLIKELY (!ml_isfinite (ML_TENSOR_AT (A, i, j)))) return ML_ERR_NAN_INPUT;
    nn = (size_t) n * (size_t) n;
    need = 6 * nn * sizeof (double);
    if (ML_UNLIKELY (ws->capacity < need)) return ML_ERR_WORKSPACE;
    buf = (double *) ml_workspace_alloc (ws, need);
    if (ML_UNLIKELY (!buf)) return ML_ERR_WORKSPACE;
    if (mat_sqrt_db (A.data, R.data, n, buf) != 0) return ML_ERR_CONVERGENCE;
    return ML_SUCCESS;
}
ML_API ml_status_t ml_matrix_exp (ml_tensor_view_t A, ml_tensor_view_t E,
                                  ml_workspace_t *ws) {
    int n = A.rows, i, j, k, s, si;
    size_t nn, need;
    double *X, *Ep, *T, *U, norm, sc;
    if (ML_UNLIKELY (n <= 0 || !A.data || !E.data || !ws)) return ML_ERR_INVALID_ARG;
    if (ML_UNLIKELY (A.cols != n || E.rows != n || E.cols != n)) return ML_ERR_INVALID_ARG;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            if (ML_UNLIKELY (!ml_isfinite (ML_TENSOR_AT (A, i, j)))) return ML_ERR_NAN_INPUT;
    nn = (size_t) n * (size_t) n;
    need = 4 * nn * sizeof (double);
    if (ML_UNLIKELY (ws->capacity < need)) return ML_ERR_WORKSPACE;
    X = (double *) ml_workspace_alloc (ws, nn * sizeof (double));
    Ep = (double *) ml_workspace_alloc (ws, nn * sizeof (double));
    T = (double *) ml_workspace_alloc (ws, nn * sizeof (double));
    U = (double *) ml_workspace_alloc (ws, nn * sizeof (double));
    if (ML_UNLIKELY (!X || !Ep || !T || !U)) return ML_ERR_WORKSPACE;
    for (i = 0; i < (int) nn; i++) X [i] = A.data [i];
    norm = mm_norm1 (X, n);
    s = 0;
    if (norm > 0.5) {
        s = (int) ml_ceil (ml_log (norm) / 0.6931471805599453) + 1;
        if (s < 0) s = 0;
        if (s > 60) s = 60;
    }
    sc = 1.0;
    for (si = 0; si < s; si++) sc *= 0.5;
    for (i = 0; i < (int) nn; i++) X [i] *= sc;
    mm_identity (Ep, n);
    for (i = 0; i < (int) nn; i++) T [i] = Ep [i];
    for (k = 1; k <= 60; k++) {
        mm_mul (U, T, X, n);
        for (i = 0; i < (int) nn; i++) T [i] = U [i] / (double) k;
        for (i = 0; i < (int) nn; i++) Ep [i] += T [i];
        if (mm_norm1 (T, n) < 1e-18 * mm_norm1 (Ep, n)) break;
    }
    for (si = 0; si < s; si++) {
        mm_mul (U, Ep, Ep, n);
        for (i = 0; i < (int) nn; i++) Ep [i] = U [i];
    }
    for (i = 0; i < (int) nn; i++) E.data [i] = Ep [i];
    return ML_SUCCESS;
}
ML_API ml_status_t ml_matrix_log (ml_tensor_view_t A, ml_tensor_view_t L,
                                  ml_workspace_t *ws) {
    int n = A.rows, i, j, k, m;
    size_t nn, need;
    double *M, *R, *S, *P, *T, *buf;
    double factor = 1.0;
    if (ML_UNLIKELY (n <= 0 || !A.data || !L.data || !ws)) return ML_ERR_INVALID_ARG;
    if (ML_UNLIKELY (A.cols != n || L.rows != n || L.cols != n)) return ML_ERR_INVALID_ARG;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            if (ML_UNLIKELY (!ml_isfinite (ML_TENSOR_AT (A, i, j)))) return ML_ERR_NAN_INPUT;
    nn = (size_t) n * (size_t) n;
    /* M (current), R (sqrt), plus 6-block db scratch, plus S,P,T,U series. */
    need = 11 * nn * sizeof (double);
    if (ML_UNLIKELY (ws->capacity < need)) return ML_ERR_WORKSPACE;
    M = (double *) ml_workspace_alloc (ws, nn * sizeof (double));
    R = (double *) ml_workspace_alloc (ws, nn * sizeof (double));
    buf = (double *) ml_workspace_alloc (ws, 6 * nn * sizeof (double));
    S = (double *) ml_workspace_alloc (ws, nn * sizeof (double));
    P = (double *) ml_workspace_alloc (ws, nn * sizeof (double));
    T = (double *) ml_workspace_alloc (ws, nn * sizeof (double));
    if (ML_UNLIKELY (!M || !R || !buf || !S || !P || !T)) return ML_ERR_WORKSPACE;
    for (i = 0; i < (int) nn; i++) M [i] = A.data [i];
    for (k = 0; k < 60; k++) {
        double dev = 0.0;
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                double d = M [(size_t) i * (size_t) n + (size_t) j] - ((i == j) ? 1.0 : 0.0);
                dev += d * d;
            }
        }
        if (ml_sqrt (dev) <= 0.5) break;
        if (mat_sqrt_db (M, R, n, buf) != 0) return ML_ERR_CONVERGENCE;
        for (i = 0; i < (int) nn; i++) M [i] = R [i];
        factor *= 2.0;
    }
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            S [(size_t) i * (size_t) n + (size_t) j] = M [(size_t) i * (size_t) n + (size_t) j] - ((i == j) ? 1.0 : 0.0);
    /* log(I+S) = S - S^2/2 + S^3/3 - ... ; result accumulates in P. */
    for (i = 0; i < (int) nn; i++) {
        P [i] = S [i];
        T [i] = S [i];
    }
    for (m = 2; m <= 100; m++) {
        double nrm = 0.0;
        mm_mul (buf, T, S, n);
        for (i = 0; i < (int) nn; i++) T [i] = buf [i];
        {
            double sgn = (m % 2 == 0) ? -1.0 : 1.0;
            for (i = 0; i < (int) nn; i++) {
                double term = sgn * T [i] / (double) m;
                P [i] += term;
                nrm += ml_fabs (term);
            }
        }
        if (nrm < 1e-18) break;
    }
    for (i = 0; i < (int) nn; i++) L.data [i] = factor * P [i];
    return ML_SUCCESS;
}
ML_API ml_status_t ml_matrix_pow (ml_tensor_view_t A, double p,
                                  ml_tensor_view_t P, ml_workspace_t *ws) {
    int n = A.rows, i;
    size_t nn, need;
    double *L;
    ml_status_t st;
    if (ML_UNLIKELY (n <= 0 || !A.data || !P.data || !ws)) return ML_ERR_INVALID_ARG;
    if (ML_UNLIKELY (A.cols != n || P.rows != n || P.cols != n)) return ML_ERR_INVALID_ARG;
    nn = (size_t) n * (size_t) n;
    need = nn * sizeof (double);
    if (ML_UNLIKELY (ws->capacity < need)) return ML_ERR_WORKSPACE;
    L = (double *) ml_workspace_alloc (ws, need);
    if (ML_UNLIKELY (!L)) return ML_ERR_WORKSPACE;
    st = ml_matrix_log (A, ml_tensor_view (L, n, n), ws);
    if (st != ML_SUCCESS) return st;
    for (i = 0; i < (int) nn; i++) L [i] *= p;
    return ml_matrix_exp (ml_tensor_view (L, n, n), P, ws);
}
