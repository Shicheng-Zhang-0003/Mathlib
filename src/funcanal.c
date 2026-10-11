#include "ml_compiler.h"
#include "ml_funcanal.h"
#include "ml_core.h"
#include "ml_spectral.h"
#include "ml_tensor.h"
#include <stddef.h>
#include <stdlib.h>

/* Banach contractions, Neumann series, Hilbert projection, Rayleigh max. */

ML_API ml_status_t ml_banach_fixed_point(ml_func_t T, double x0, double L,
                                         double tol, int max_iter,
                                         double *xstar, double *bound_out) {
    double x = x0;
    double x1;
    double first_step;
    int n = 0;
    if (ML_UNLIKELY(!T || !xstar)) return ML_ERR_INVALID_ARG;
    if (ml_isnan(x0) || !ml_isfinite(x0)) return ML_ERR_NAN_INPUT;
    if (ml_isnan(L) || !ml_isfinite(L) || L < 0.0 || L >= 1.0) return ML_ERR_INVALID_ARG;
    if (!(tol > 0.0) || !ml_isfinite(tol)) tol = 1e-12;
    if (max_iter <= 0) max_iter = 10000;
    if (max_iter > 1000000) max_iter = 1000000;
    x1 = T(x0);
    if (ml_isnan(x1) || !ml_isfinite(x1)) return ML_ERR_NAN_INPUT;
    first_step = ml_fabs(x1 - x0);
    x = x1;
    for (n = 1; n <= max_iter; n++) {
        double xn = T(x);
        if (ml_isnan(xn) || !ml_isfinite(xn)) return ML_ERR_NAN_INPUT;
        /* A posteriori bound from the certificate, reported at exit. */
        if (ml_fabs(xn - x) <= tol * (1.0 + ml_fabs(xn))) {
            long double pw = 1.0L;
            int i;
            for (i = 0; i < n; i++) pw *= (long double)L;
            *xstar = xn;
            if (bound_out) {
                long double b = (L >= 1.0L - 1e-300L) ? 1e300L
                    : pw / (1.0L - (long double)L) * (long double)first_step;
                double bo = (double)b;
                *bound_out = ml_isfinite(bo) ? bo : ml_make_nan();
            }
            return ML_SUCCESS;
        }
        x = xn;
    }
    *xstar = ml_make_nan();
    if (bound_out) *bound_out = ml_make_nan();
    return ML_ERR_CONVERGENCE;
}

ML_API ml_status_t ml_neumann_inverse(const double *A, int n, int terms, double *out) {
    /* S = I + B + B^2 + ... with B = I - A; needs ||B||_inf < 1. n <= 16. */
    double b[256];
    double acc[256];
    double pw[256];
    int i;
    int j;
    int t;
    if (!A || !out || n <= 0 || n > 16) return ML_ERR_INVALID_ARG;
    if (terms <= 0 || terms > 1000) return ML_ERR_INVALID_ARG;
    for (i = 0; i < n; i++) {
        long double rs = 0.0L;
        for (j = 0; j < n; j++) {
            double a = A[(size_t)i * (size_t)n + (size_t)j];
            double bij = ((i == j) ? 1.0 : 0.0) - a;
            if (ml_isnan(a) || !ml_isfinite(a)) return ML_ERR_NAN_INPUT;
            b[(size_t)i * (size_t)n + (size_t)j] = bij;
            rs += (bij < 0.0 ? -bij : bij);
        }
        if (!(rs < 1.0L)) return ML_ERR_INVALID_ARG;
    }
    for (i = 0; i < n * n; i++) {
        acc[i] = (i / n == i % n) ? 1.0 : 0.0;
        pw[i] = acc[i];
    }
    for (t = 1; t < terms; t++) {
        double nxt[256];
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                long double s = 0.0L;
                int k;
                for (k = 0; k < n; k++) s += (long double)pw[(size_t)i * (size_t)n + (size_t)k] * (long double)b[(size_t)k * (size_t)n + (size_t)j];
                nxt[(size_t)i * (size_t)n + (size_t)j] = (double)s;
                if (!ml_isfinite(nxt[(size_t)i * (size_t)n + (size_t)j])) return ML_ERR_OVERFLOW;
            }
        }
        for (i = 0; i < n * n; i++) {
            acc[i] += nxt[i];
            if (!ml_isfinite(acc[i])) return ML_ERR_OVERFLOW;
            pw[i] = nxt[i];
        }
    }
    for (i = 0; i < n * n; i++) out[i] = acc[i];
    return ML_SUCCESS;
}

ML_API ml_status_t ml_hilbert_project(const double *V, const double *y,
                                      int m, int k, double *coeff) {
    /* Modified Gram-Schmidt on columns, then back substitution on R c = Q'y.
     * Rank loss (tiny R_ii) is SINGULAR, never silent. */
    double q[32768];
    double r[1024];
    double qty[32];
    int i;
    int j;
    int t;
    if (!V || !y || !coeff || m <= 0 || m > 1024 || k <= 0 || k > 32) return ML_ERR_INVALID_ARG;
    for (i = 0; i < m; i++) {
        if (ml_isnan(y[i]) || !ml_isfinite(y[i])) return ML_ERR_NAN_INPUT;
        for (j = 0; j < k; j++) {
            double v = V[(size_t)i * (size_t)k + (size_t)j];
            if (ml_isnan(v) || !ml_isfinite(v)) return ML_ERR_NAN_INPUT;
            q[(size_t)i * (size_t)k + (size_t)j] = v;
        }
    }
    for (j = 0; j < k; j++) r[j] = 0.0;
    for (j = 0; j < k; j++) {
        long double colnrm = 0.0L;
        for (t = 0; t < m; t++) {
            long double v = (long double)q[(size_t)t * (size_t)k + (size_t)j];
            colnrm += v * v;
        }
        colnrm = __builtin_sqrtl(colnrm);
        for (i = 0; i < j; i++) {
            long double s = 0.0L;
            for (t = 0; t < m; t++) s += (long double)q[(size_t)t * (size_t)k + (size_t)j] * (long double)q[(size_t)t * (size_t)k + (size_t)i];
            r[(size_t)i * (size_t)k + (size_t)j] = (double)s;
            for (t = 0; t < m; t++) q[(size_t)t * (size_t)k + (size_t)j] -= (double)s * q[(size_t)t * (size_t)k + (size_t)i];
        }
        {
            long double nrm = 0.0L;
            for (t = 0; t < m; t++) {
                long double v = (long double)q[(size_t)t * (size_t)k + (size_t)j];
                nrm += v * v;
            }
            nrm = __builtin_sqrtl(nrm);
            /* Scale-aware rank test: irrational norms never hit exact zero
             * (collinear columns leave ~1e-19 residue, not 0). */
            if (!(nrm > 1e-12L * (1.0L + colnrm)) || !ml_isfinite((double)nrm))
                return ML_ERR_SINGULAR;
            r[(size_t)j * (size_t)k + (size_t)j] = (double)nrm;
            for (t = 0; t < m; t++) q[(size_t)t * (size_t)k + (size_t)j] /= (double)nrm;
        }
    }
    for (i = 0; i < k; i++) {
        long double s = 0.0L;
        for (t = 0; t < m; t++) s += (long double)q[(size_t)t * (size_t)k + (size_t)i] * (long double)y[t];
        qty[i] = (double)s;
    }
    for (i = k - 1; i >= 0; i--) {
        long double s = (long double)qty[i];
        int j2;
        double d = r[(size_t)i * (size_t)k + (size_t)i];
        for (j2 = i + 1; j2 < k; j2++) s -= (long double)r[(size_t)i * (size_t)k + (size_t)j2] * (long double)coeff[j2];
        if (d == 0.0) return ML_ERR_SINGULAR;
        coeff[i] = (double)(s / (long double)d);
        if (!ml_isfinite(coeff[i])) return ML_ERR_OVERFLOW;
    }
    return ML_SUCCESS;
}

ML_API ml_status_t ml_rayleigh_max(const double *A, int n, double *lambda) {
    /* Thin wrapper over the staged power iteration (row-major in). */
    double *buf = NULL;
    ml_tensor_view_t view;
    double lam = 0.0;
    ml_status_t st;
    double vec[128];
    if (!A || !lambda || n <= 0 || n > 128) return ML_ERR_INVALID_ARG;
    {
        int i;
        for (i = 0; i < n * n; i++) {
            if (ml_isnan(A[i]) || !ml_isfinite(A[i])) return ML_ERR_NAN_INPUT;
        }
    }
    if ((size_t)n * (size_t)n > (size_t)-1 / sizeof(double)) return ML_ERR_INVALID_ARG;
    buf = (double *)malloc((size_t)n * (size_t)n * sizeof(double));
    if (!buf) return ML_ERR_WORKSPACE;
    {
        int i;
        for (i = 0; i < n * n; i++) buf[i] = A[i];
    }
    view.data = buf;
    view.rows = n;
    view.cols = n;
    st = ml_power_iter(view, &lam, vec, 0, 0.0);
    free(buf);
    if (st != ML_SUCCESS) return st;
    *lambda = lam;
    return ML_SUCCESS;
}
