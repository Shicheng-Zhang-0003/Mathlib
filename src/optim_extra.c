#include "ml_compiler.h"
#include "ml_optim_extra.h"

/* ============================================================================
 * v12R2 optimisation depth pack.
 * ========================================================================== */

#define OE_EPS 1e-9
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
                double f = A [(size_t) i * (size_t) n + (size_t) k] / d;
                b [i] -= f * b [k];
                for (j = k; j < n; j++)
                    A [(size_t) i * (size_t) n + (size_t) j] -=
                        f * A [(size_t) k * (size_t) n + (size_t) j];
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
/* --------------------------------------------------------- simplex ------- */

#define OE_M 32
#define OE_N 32
static ml_status_t oe_simplex (const double *c, const double *A, const double *b,
                               int m, int n, double *x, double *obj) {
    double tab [(OE_M + 1) * (OE_N + OE_M + 1)];
    int basis [OE_M];
    int cols = n + m + 1, i, j, iter;
    if (m <= 0 || n <= 0 || m > OE_M || n > OE_N) return ML_ERR_INVALID_ARG;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++)
            tab [(size_t) i * (size_t) cols + (size_t) j] = A [(size_t) i * (size_t) n + (size_t) j];
        for (j = 0; j < m; j++)
            tab [(size_t) i * (size_t) cols + (size_t) (n + j)] = (i == j) ? 1.0 : 0.0;
        tab [(size_t) i * (size_t) cols + (size_t) (cols - 1)] = b [i];
        basis [i] = n + i;
        if (b [i] < 0.0) return ML_ERR_INVALID_ARG;
    }
    for (j = 0; j < cols; j++) tab [(size_t) m * (size_t) cols + (size_t) j] = 0.0;
    for (j = 0; j < n; j++)
        tab [(size_t) m * (size_t) cols + (size_t) j] = -c [j];
    for (iter = 0; iter < 10000; iter++) {
        int enter = -1, leave = -1;
        double bestratio = 0.0;
        /* Bland: smallest index with negative reduced cost. */
        for (j = 0; j < cols - 1; j++)
            if (tab [(size_t) m * (size_t) cols + (size_t) j] < -OE_EPS) {
                enter = j;
                break;
            }
        if (enter < 0) break;
        for (i = 0; i < m; i++) {
            double a = tab [(size_t) i * (size_t) cols + (size_t) enter];
            if (a > OE_EPS) {
                double r = tab [(size_t) i * (size_t) cols + (size_t) (cols - 1)] / a;
                if (leave < 0 || r < bestratio - OE_EPS ||
                    (r < bestratio + OE_EPS && i < leave)) {
                    bestratio = r;
                    leave = i;
                }
            }
        }
        if (leave < 0) return ML_ERR_CONVERGENCE; /* unbounded */
        {
            double piv = tab [(size_t) leave * (size_t) cols + (size_t) enter];
            basis [leave] = enter;
            for (j = 0; j < cols; j++)
                tab [(size_t) leave * (size_t) cols + (size_t) j] /= piv;
            for (i = 0; i < m + 1; i++) {
                double f;
                if (i == leave) continue;
                f = tab [(size_t) i * (size_t) cols + (size_t) enter];
                if (f == 0.0) continue;
                for (j = 0; j < cols; j++)
                    tab [(size_t) i * (size_t) cols + (size_t) j] -=
                        f * tab [(size_t) leave * (size_t) cols + (size_t) j];
            }
        }
    }
    for (j = 0; j < n; j++) x [j] = 0.0;
    for (i = 0; i < m; i++)
        if (basis [i] >= 0 && basis [i] < n)
            x [basis [i]] = tab [(size_t) i * (size_t) cols + (size_t) (cols - 1)];
    if (obj) *obj = tab [(size_t) m * (size_t) cols + (size_t) (cols - 1)];
    return ML_SUCCESS;
}
ML_API ml_status_t ml_simplex_lp (const double *c, const double *A,
                                  const double *b, int m, int n, double *x,
                                  double *obj) {
    if (ML_UNLIKELY (!c || !A || !b || !x)) return ML_ERR_INVALID_ARG;
    return oe_simplex (c, A, b, m, n, x, obj);
}
ML_API ml_status_t ml_simplex_lp_min (const double *c, const double *A,
                                      const double *b, int m, int n, double *x,
                                      double *obj) {
    double nc [OE_N];
    int j, mr = 0;
    if (ML_UNLIKELY (!c || !A || !b || !x || n > OE_N)) return ML_ERR_INVALID_ARG;
    for (j = 0; j < n; j++) nc [j] = -c [j];
    mr = oe_simplex (nc, A, b, m, n, x, obj);
    if (mr == ML_SUCCESS && obj) *obj = -*obj;
    return mr;
}
/* ---------------------------------------------- branch & bound ILP ------- */

static void oe_bb (const double *c, double *A, double *b, int m, int n,
                   double *x, double *best, double *bestobj, int *nodes,
                   int max_nodes) {
    double xr [OE_N], objr = 0.0;
    int i, j, frac = -1;
    if (++(*nodes) > max_nodes) return;
    if (oe_simplex (c, A, b, m, n, xr, &objr) != ML_SUCCESS) return;
    if (objr <= *bestobj + 1e-9) return;
    for (j = 0; j < n; j++) {
        double d = xr [j] - (double) (long) (xr [j] + 0.5);
        if (ml_fabs (d) > 1e-6) {
            frac = j;
            break;
        }
    }
    if (frac < 0) {
        *bestobj = objr;
        for (j = 0; j < n; j++) best [j] = xr [j];
        return;
    }
    if (m + 2 > OE_M) return;
    {
        double fl = ml_floor (xr [frac]);
        /* branch x_frac <= fl */
        for (j = 0; j < n; j++) A [(size_t) m * (size_t) n + (size_t) j] = 0.0;
        A [(size_t) m * (size_t) n + (size_t) frac] = 1.0;
        b [m] = fl;
        oe_bb (c, A, b, m + 1, n, x, best, bestobj, nodes, max_nodes);
        /* branch x_frac >= fl+1  =>  -x_frac <= -(fl+1) */
        for (j = 0; j < n; j++) A [(size_t) m * (size_t) n + (size_t) j] = 0.0;
        A [(size_t) m * (size_t) n + (size_t) frac] = -1.0;
        b [m] = -(fl + 1.0);
        oe_bb (c, A, b, m + 1, n, x, best, bestobj, nodes, max_nodes);
    }
    (void) i;
}
ML_API ml_status_t ml_integer_lp (const double *c, const double *A,
                                  const double *b, int m, int n, double *x,
                                  double *obj, int max_nodes) {
    double Ab [OE_M * OE_N], bb [OE_M], best [OE_N], bestobj = -HUGE_VAL;
    int i, j, nodes = 0;
    if (ML_UNLIKELY (!c || !A || !b || !x || m < 0 || n <= 0 || m > OE_M - 2))
        return ML_ERR_INVALID_ARG;
    for (i = 0; i < m; i++) {
        bb [i] = b [i];
        for (j = 0; j < n; j++)
            Ab [(size_t) i * (size_t) n + (size_t) j] = A [(size_t) i * (size_t) n + (size_t) j];
    }
    if (max_nodes <= 0) max_nodes = 10000;
    for (j = 0; j < n; j++) best [j] = 0.0;
    oe_bb (c, Ab, bb, m, n, x, best, &bestobj, &nodes, max_nodes);
    if (bestobj == -HUGE_VAL) return ML_ERR_SINGULAR;
    for (j = 0; j < n; j++) x [j] = best [j];
    if (obj) *obj = bestobj;
    return ML_SUCCESS;
}
/* ------------------------------------------------- quadratic prog -------- */

ML_API ml_status_t ml_qp_eq (const double *H, const double *g, const double *A,
                             const double *b, int n, int m, double *x) {
    double K [64 * 64], rhs [64], sol [64];
    int i, j, N;
    if (ML_UNLIKELY (!H || !g || !x || n <= 0 || m < 0 || n + m > 64))
        return ML_ERR_INVALID_ARG;
    if (m > 0 && (!A || !b)) return ML_ERR_INVALID_ARG;
    N = n + m;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            if (i < n && j < n) K [(size_t) i * (size_t) N + (size_t) j] = H [(size_t) i * (size_t) n + (size_t) j];
            else if (i < n && j >= n) K [(size_t) i * (size_t) N + (size_t) j] = A [(size_t) (j - n) * (size_t) n + (size_t) i];
            else if (i >= n && j < n) K [(size_t) i * (size_t) N + (size_t) j] = A [(size_t) (i - n) * (size_t) n + (size_t) j];
            else K [(size_t) i * (size_t) N + (size_t) j] = 0.0;
        }
        rhs [i] = (i < n) ? -g [i] : b [i - n];
    }
    if (oe_solve (K, rhs, sol, N) != ML_SUCCESS) return ML_ERR_SINGULAR;
    for (i = 0; i < n; i++) x [i] = sol [i];
    return ML_SUCCESS;
}
static ml_status_t oe_qp_eq_ws (const double *H, const double *g, const double *x,
                                const int *ws, int nws, const double *A, int n,
                                double *p, double *lam) {
    double K [64 * 64], rhs [64], sol [64];
    double Hx [OE_N];
    int i, j, N = n + nws;
    if (N > 64) return ML_ERR_INVALID_ARG;
    for (i = 0; i < n; i++) {
        double s = 0.0;
        for (j = 0; j < n; j++) s += H [(size_t) i * (size_t) n + (size_t) j] * x [j];
        Hx [i] = s + g [i];
    }
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            if (i < n && j < n) K [(size_t) i * (size_t) N + (size_t) j] = H [(size_t) i * (size_t) n + (size_t) j];
            else if (i < n) K [(size_t) i * (size_t) N + (size_t) j] = A [(size_t) ws [j - n] * (size_t) n + (size_t) i];
            else if (j < n) K [(size_t) i * (size_t) N + (size_t) j] = A [(size_t) ws [i - n] * (size_t) n + (size_t) j];
            else K [(size_t) i * (size_t) N + (size_t) j] = 0.0;
        }
        rhs [i] = (i < n) ? -Hx [i] : 0.0;
    }
    if (oe_solve (K, rhs, sol, N) != ML_SUCCESS) return ML_ERR_SINGULAR;
    for (i = 0; i < n; i++) p [i] = sol [i];
    for (i = 0; i < nws; i++) lam [i] = sol [n + i];
    return ML_SUCCESS;
}
ML_API ml_status_t ml_qp_active_set (const double *H, const double *g,
                                     const double *A, const double *b, int n,
                                     int m, const double *x0, int max_iter,
                                     double tol, double *x) {
    int ws [OE_M], nws = 0, it, i, j;
    double lam [OE_M];
    if (ML_UNLIKELY (!H || !g || !A || !b || !x || n <= 0 || n > OE_N || m < 0 || m > OE_M))
        return ML_ERR_INVALID_ARG;
    if (max_iter <= 0) max_iter = 200;
    if (tol <= 0.0) tol = 1e-9;
    for (j = 0; j < n; j++) x [j] = x0 ? x0 [j] : 0.0;
    for (it = 0; it < max_iter; it++) {
        double p [OE_N], pn = 0.0;
        int blk;
        if (oe_qp_eq_ws (H, g, x, ws, nws, A, n, p, lam) != ML_SUCCESS)
            return ML_ERR_SINGULAR;
        for (j = 0; j < n; j++) pn += p [j] * p [j];
        pn = ml_sqrt (pn);
        if (pn <= tol) {
            int minlam = -1;
            double mlam = 0.0;
            for (i = 0; i < nws; i++)
                if (lam [i] < mlam) {
                    mlam = lam [i];
                    minlam = i;
                }
            if (minlam < 0) {
                /* feasible and optimal */
                return ML_SUCCESS;
            }
            for (i = minlam; i < nws - 1; i++) ws [i] = ws [i + 1];
            nws--;
            continue;
        }
        /* find max step to a violated constraint */
        blk = -1;
        {
            double alpha = 1.0;
            for (i = 0; i < m; i++) {
                int inws = 0;
                double ap = 0.0, slack;
                for (j = 0; j < nws; j++)
                    if (ws [j] == i) inws = 1;
                if (inws) continue;
                for (j = 0; j < n; j++) ap += A [(size_t) i * (size_t) n + (size_t) j] * p [j];
                if (ap <= OE_EPS) continue;
                slack = b [i];
                for (j = 0; j < n; j++) slack -= A [(size_t) i * (size_t) n + (size_t) j] * x [j];
                if (slack / ap < alpha) {
                    alpha = slack / ap;
                    blk = i;
                }
            }
            if (alpha >= 1.0) {
                for (j = 0; j < n; j++) x [j] += p [j];
            } else {
                for (j = 0; j < n; j++) x [j] += alpha * p [j];
                if (nws >= OE_M) return ML_ERR_CONVERGENCE;
                ws [nws++] = blk;
            }
        }
    }
    return ML_ERR_CONVERGENCE;
}
/* ------------------------------------------------ internal BFGS ---------- */

typedef double (*oe_fn_t) (const double *x, int n, void *ctx);
static double oe_grad (oe_fn_t f, void *ctx, const double *x, int n, double *g) {
    int i;
    double fx = f (x, n, ctx);
    for (i = 0; i < n; i++) {
        double xi = x [i], h = 1e-7 * (ml_fabs (xi) + 1.0);
        double xp [OE_N], xm [OE_N];
        int k;
        for (k = 0; k < n; k++) {
            xp [k] = x [k];
            xm [k] = x [k];
        }
        xp [i] = xi + h;
        xm [i] = xi - h;
        g [i] = (f (xp, n, ctx) - f (xm, n, ctx)) / (2.0 * h);
    }
    return fx;
}
static ml_status_t oe_bfgs (oe_fn_t f, void *ctx, double *x, int n,
                            int max_iter, double tol, double *fout) {
    double Hinv [OE_N * OE_N], g [OE_N], gn [OE_N], p [OE_N], s [OE_N], y [OE_N];
    int i, j, it;
    if (n <= 0 || n > OE_N) return ML_ERR_INVALID_ARG;
    if (max_iter <= 0) max_iter = 500;
    if (tol <= 0.0) tol = 1e-10;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++) Hinv [(size_t) i * (size_t) n + (size_t) j] = (i == j) ? 1.0 : 0.0;
    oe_grad (f, ctx, x, n, g);
    for (it = 0; it < max_iter; it++) {
        double gnorm = 0.0, alpha, fy;
        for (j = 0; j < n; j++) gnorm += g [j] * g [j];
        gnorm = ml_sqrt (gnorm);
        if (gnorm <= tol) break;
        for (j = 0; j < n; j++) {
            double s2 = 0.0;
            for (i = 0; i < n; i++) s2 -= Hinv [(size_t) j * (size_t) n + (size_t) i] * g [i];
            p [j] = s2;
        }
        /* Armijo backtracking */
        alpha = 1.0;
        {
            double gp = 0.0, f0 = f (x, n, ctx);
            for (j = 0; j < n; j++) gp += g [j] * p [j];
            for (;;) {
                double xt [OE_N];
                for (j = 0; j < n; j++) xt [j] = x [j] + alpha * p [j];
                fy = f (xt, n, ctx);
                if (fy <= f0 + 1e-4 * alpha * gp || alpha < 1e-12) break;
                alpha *= 0.5;
            }
        }
        for (j = 0; j < n; j++) {
            s [j] = alpha * p [j];
            x [j] += s [j];
        }
        {
            double ys = 0.0;
            oe_grad (f, ctx, x, n, gn);
            for (j = 0; j < n; j++) {
                y [j] = gn [j] - g [j];
                ys += y [j] * s [j];
            }
            if (ys > 1e-14) {
                double rho = 1.0 / ys;
                double Hp [OE_N * OE_N];
                /* Hinv <- (I - rho s y^T) Hinv (I - rho y s^T) + rho s s^T */
                for (i = 0; i < n; i++)
                    for (j = 0; j < n; j++) {
                        double v = 0.0;
                        int a, b;
                        for (a = 0; a < n; a++)
                            for (b = 0; b < n; b++) {
                                double lia = (i == a) ? 1.0 : 0.0;
                                double rib = (j == b) ? 1.0 : 0.0;
                                double m1 = lia - rho * s [i] * y [a];
                                double m2 = rib - rho * y [b] * s [j];
                                v += m1 * Hinv [(size_t) a * (size_t) n + (size_t) b] * m2;
                            }
                        Hp [(size_t) i * (size_t) n + (size_t) j] = v + rho * s [i] * s [j];
                    }
                for (i = 0; i < n * n; i++) Hinv [i] = Hp [i];
            }
            for (j = 0; j < n; j++) g [j] = gn [j];
        }
    }
    if (fout) *fout = f (x, n, ctx);
    return ML_SUCCESS;
}
/* ------------------------------------------------- nonlinear LS ---------- */

ML_API ml_status_t ml_gauss_newton (ml_residual_t res, int m, const double *x0,
                                    int n, int max_iter, double tol,
                                    double *x_out) {
    double J [(128) * OE_N], JTJ [OE_N * OE_N], JTr [OE_N], step [OE_N], r [128], rp [128];
    int i, j, k, it;
    if (ML_UNLIKELY (!res || !x0 || !x_out || n <= 0 || n > OE_N || m <= 0 || m > 128))
        return ML_ERR_INVALID_ARG;
    if (max_iter <= 0) max_iter = 100;
    if (tol <= 0.0) tol = 1e-10;
    for (j = 0; j < n; j++) x_out [j] = x0 [j];
    for (it = 0; it < max_iter; it++) {
        double ss = 0.0, gnorm;
        res (x_out, n, r, NULL);
        for (i = 0; i < m; i++) ss += r [i] * r [i];
        /* finite-difference Jacobian */
        for (j = 0; j < n; j++) {
            double xj = x_out [j], h = 1e-7 * (ml_fabs (xj) + 1.0);
            x_out [j] = xj + h;
            res (x_out, n, rp, NULL);
            x_out [j] = xj;
            for (i = 0; i < m; i++)
                J [(size_t) i * (size_t) n + (size_t) j] = (rp [i] - r [i]) / h;
        }
        for (i = 0; i < n; i++) {
            JTr [i] = 0.0;
            for (k = 0; k < m; k++) JTr [i] += J [(size_t) k * (size_t) n + (size_t) i] * r [k];
            for (j = 0; j < n; j++) {
                double s = 0.0;
                for (k = 0; k < m; k++)
                    s += J [(size_t) k * (size_t) n + (size_t) i] * J [(size_t) k * (size_t) n + (size_t) j];
                JTJ [(size_t) i * (size_t) n + (size_t) j] = s;
            }
        }
        gnorm = 0.0;
        for (i = 0; i < n; i++) gnorm += JTr [i] * JTr [i];
        gnorm = ml_sqrt (gnorm);
        if (gnorm <= tol) break;
        for (i = 0; i < n; i++) {
            JTJ [(size_t) i * (size_t) n + (size_t) i] *= 1.0 + 1e-9;
            step [i] = -JTr [i];
        }
        if (oe_solve (JTJ, step, step, n) != ML_SUCCESS) return ML_ERR_SINGULAR;
        {
            double alpha = 1.0, ssnew;
            for (;;) {
                double xt [OE_N], s2 = 0.0;
                int bad = 0;
                for (j = 0; j < n; j++) xt [j] = x_out [j] + alpha * step [j];
                res (xt, n, rp, NULL);
                for (i = 0; i < m; i++) s2 += rp [i] * rp [i];
                ssnew = s2;
                if (ssnew <= ss || alpha < 1e-10) {
                    (void) bad;
                    break;
                }
                alpha *= 0.5;
            }
            for (j = 0; j < n; j++) x_out [j] += alpha * step [j];
        }
    }
    return ML_SUCCESS;
}
/* ------------------------------------------------- Newton / boxes -------- */

ML_API ml_status_t ml_newton_min (ml_vec_func_t f, const double *x0, int n,
                                  int max_iter, double tol, double *x_out,
                                  double *f_out) {
    double x [OE_N], g [OE_N], Hs [OE_N * OE_N], step [OE_N];
    int i, j, k, it;
    if (ML_UNLIKELY (!f || !x0 || !x_out || n <= 0 || n > OE_N))
        return ML_ERR_INVALID_ARG;
    if (max_iter <= 0) max_iter = 200;
    if (tol <= 0.0) tol = 1e-10;
    for (j = 0; j < n; j++) x [j] = x0 [j];
    for (it = 0; it < max_iter; it++) {
        double gnorm = 0.0;
        for (i = 0; i < n; i++) {
            double xi = x [i], h = 1e-6 * (ml_fabs (xi) + 1.0);
            double xp [OE_N], xm [OE_N];
            for (k = 0; k < n; k++) {
                xp [k] = x [k];
                xm [k] = x [k];
            }
            xp [i] = xi + h;
            xm [i] = xi - h;
            g [i] = (f (xp, n) - f (xm, n)) / (2.0 * h);
        }
        for (i = 0; i < n; i++) {
            double xi = x [i], hi = 1e-4 * (ml_fabs (xi) + 1.0);
            for (j = 0; j < n; j++) {
                if (i == j) {
                    double xp [OE_N], xm [OE_N], f0 = f (x, n);
                    for (k = 0; k < n; k++) {
                        xp [k] = x [k];
                        xm [k] = x [k];
                    }
                    xp [i] = xi + hi;
                    xm [i] = xi - hi;
                    Hs [(size_t) i * (size_t) n + (size_t) j] =
                        (f (xp, n) - 2.0 * f0 + f (xm, n)) / (hi * hi);
                } else {
                    double xpp [OE_N], xpm [OE_N], xmp [OE_N], xmm [OE_N];
                    double hj = 1e-4 * (ml_fabs (x [j]) + 1.0);
                    for (k = 0; k < n; k++) {
                        xpp [k] = x [k];
                        xpm [k] = x [k];
                        xmp [k] = x [k];
                        xmm [k] = x [k];
                    }
                    xpp [i] += hi;
                    xpp [j] += hj;
                    xpm [i] += hi;
                    xpm [j] -= hj;
                    xmp [i] -= hi;
                    xmp [j] += hj;
                    xmm [i] -= hi;
                    xmm [j] -= hj;
                    Hs [(size_t) i * (size_t) n + (size_t) j] =
                        (f (xpp, n) - f (xpm, n) - f (xmp, n) + f (xmm, n)) / (4.0 * hi * hj);
                }
            }
        }
        for (j = 0; j < n; j++) {
            gnorm += g [j] * g [j];
            step [j] = -g [j];
        }
        gnorm = ml_sqrt (gnorm);
        if (gnorm <= tol) break;
        {
            double K [OE_N * OE_N];
            for (i = 0; i < n * n; i++) K [i] = Hs [i];
            for (i = 0; i < n; i++) K [(size_t) i * (size_t) n + (size_t) i] += 1e-8;
            if (oe_solve (K, step, step, n) != ML_SUCCESS) return ML_ERR_SINGULAR;
        }
        {
            double f0 = f (x, n), alpha = 1.0, gp = 0.0;
            for (j = 0; j < n; j++) gp += g [j] * step [j];
            for (;;) {
                double xt [OE_N];
                for (j = 0; j < n; j++) xt [j] = x [j] + alpha * step [j];
                if (f (xt, n) <= f0 + 1e-4 * alpha * gp || alpha < 1e-12) break;
                alpha *= 0.5;
            }
            for (j = 0; j < n; j++) x [j] += alpha * step [j];
        }
    }
    for (j = 0; j < n; j++) x_out [j] = x [j];
    if (f_out) *f_out = f (x, n);
    return ML_SUCCESS;
}
static double oe_clamp (double v, double lo, double hi) {
    if (v < lo) return lo;
    if (v > hi) return hi;
    return v;
}
ML_API ml_status_t ml_projected_gradient_box (ml_vec_func_t f,
                                              const double *x0, const double *lo,
                                              const double *up, int n,
                                              int max_iter, double tol,
                                              double *x_out, double *f_out) {
    double x [OE_N], g [OE_N];
    int i, j, it;
    if (ML_UNLIKELY (!f || !x0 || !lo || !up || !x_out || n <= 0 || n > OE_N))
        return ML_ERR_INVALID_ARG;
    if (max_iter <= 0) max_iter = 500;
    if (tol <= 0.0) tol = 1e-10;
    for (j = 0; j < n; j++) x [j] = oe_clamp (x0 [j], lo [j], up [j]);
    for (it = 0; it < max_iter; it++) {
        int k;
        double step = 1.0, gnorm = 0.0, fx = f (x, n);
        for (i = 0; i < n; i++) {
            double xi = x [i], h = 1e-7 * (ml_fabs (xi) + 1.0);
            double xp [OE_N], xm [OE_N];
            for (k = 0; k < n; k++) {
                xp [k] = x [k];
                xm [k] = x [k];
            }
            xp [i] = oe_clamp (xi + h, lo [i], up [i]);
            xm [i] = oe_clamp (xi - h, lo [i], up [i]);
            g [i] = (f (xp, n) - f (xm, n)) / (xp [i] - xm [i] + 1e-300);
            gnorm += g [i] * g [i];
        }
        gnorm = ml_sqrt (gnorm);
        if (gnorm <= tol) break;
        for (;;) {
            double xt [OE_N], ft;
            for (j = 0; j < n; j++)
                xt [j] = oe_clamp (x [j] - step * g [j], lo [j], up [j]);
            ft = f (xt, n);
            if (ft < fx || step < 1e-12) {
                for (j = 0; j < n; j++) x [j] = xt [j];
                break;
            }
            step *= 0.5;
        }
    }
    for (j = 0; j < n; j++) x_out [j] = x [j];
    if (f_out) *f_out = f (x, n);
    return ML_SUCCESS;
}
/* ------------------------------------------------ augmented Lagrangian --- */

typedef struct {
    ml_vec_func_t f;
    ml_residual_t h, g;
    int nh, ng;
    void *hctx, *gctx;
    double mu;
    double lamh [32];
    double lamg [32];
} oe_al_ctx_t;
static double oe_al_fn (const double *x, int n, void *vctx) {
    oe_al_ctx_t *c = (oe_al_ctx_t *) vctx;
    double h [32], g [32], val = c->f (x, n);
    int i;
    if (c->nh > 0) {
        c->h (x, n, h, c->hctx);
        for (i = 0; i < c->nh; i++) {
            double t = h [i] + c->lamh [i] / c->mu;
            val += c->mu * 0.5 * t * t - c->lamh [i] * c->lamh [i] / (2.0 * c->mu);
        }
    }
    if (c->ng > 0) {
        c->g (x, n, g, c->gctx);
        for (i = 0; i < c->ng; i++) {
            double t = c->lamg [i] + c->mu * g [i];
            if (t < 0.0) t = 0.0;
            val += (t * t - c->lamg [i] * c->lamg [i]) / (2.0 * c->mu);
        }
    }
    return val;
}
ML_API ml_status_t ml_augmented_lagrangian (ml_vec_func_t f, ml_residual_t h,
                                            int nh, ml_residual_t g, int ng,
                                            void *hctx, void *gctx,
                                            const double *x0, int n,
                                            int max_iter, double tol,
                                            double *x_out, double *f_out) {
    oe_al_ctx_t ctx;
    double x [OE_N];
    int i, it, j;
    if (ML_UNLIKELY (!f || !x0 || !x_out || n <= 0 || n > OE_N || nh < 0 || ng < 0 || nh > 32 || ng > 32))
        return ML_ERR_INVALID_ARG;
    if ((nh > 0 && !h) || (ng > 0 && !g)) return ML_ERR_INVALID_ARG;
    if (max_iter <= 0) max_iter = 30;
    if (tol <= 0.0) tol = 1e-8;
    ctx.f = f;
    ctx.h = h;
    ctx.g = g;
    ctx.nh = nh;
    ctx.ng = ng;
    ctx.hctx = hctx;
    ctx.gctx = gctx;
    ctx.mu = 1.0;
    for (i = 0; i < nh; i++) ctx.lamh [i] = 0.0;
    for (i = 0; i < ng; i++) ctx.lamg [i] = 0.0;
    for (j = 0; j < n; j++) x [j] = x0 [j];
    for (it = 0; it < max_iter; it++) {
        double cnorm = 0.0;
        oe_bfgs (oe_al_fn, &ctx, x, n, 300, 1e-8, NULL);
        if (nh > 0) {
            double hv [32];
            ctx.h (x, n, hv, hctx);
            for (i = 0; i < nh; i++) {
                ctx.lamh [i] += ctx.mu * hv [i];
                cnorm += hv [i] * hv [i];
            }
        }
        if (ng > 0) {
            double gv [32];
            ctx.g (x, n, gv, gctx);
            for (i = 0; i < ng; i++) {
                double lam = ctx.lamg [i] + ctx.mu * gv [i];
                ctx.lamg [i] = (lam > 0.0) ? lam : 0.0;
                if (gv [i] > 0.0) cnorm += gv [i] * gv [i];
            }
        }
        if (ml_sqrt (cnorm) <= tol) break;
        ctx.mu *= 10.0;
    }
    for (j = 0; j < n; j++) x_out [j] = x [j];
    if (f_out) *f_out = f (x, n);
    return ML_SUCCESS;
}
