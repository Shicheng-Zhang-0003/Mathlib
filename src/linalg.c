#include "ml_compiler.h"
#include "ml_linalg.h"
#include "internal/hypot.h"
#include <stdint.h>

/* v11S CLOSURE IP-10: linear algebra edge hardening */

ML_API ml_status_t ml_lu_decomp(ml_tensor_view_t A, ml_tensor_view_t LU, int* P, ml_workspace_t* ws) {
    (void)ws;

    int n = A.rows;

    /* SAFETY BY DEFAULT: Unconditional NULL and dimension checks */
    if (ML_UNLIKELY(n <= 0 || A.data == NULL || LU.data == NULL || P == NULL)) {
        return ML_ERR_INVALID_ARG;
    }

    if (ML_UNLIKELY(A.cols != n || LU.rows != n || LU.cols != n)) {
        return ML_ERR_INVALID_ARG;
    }

    /* size_t overflow guard for (r*cols+c) indexing via ML_TENSOR_AT.
     * On 32-bit, n~46341+ can wrap; reject before any access. */
    {
        size_t sn = (size_t)n;
        const size_t size_max = (size_t)-1;
        if (ML_UNLIKELY(sn > 0 && sn - 1 > (size_max - 1) / sn)) {
            return ML_ERR_INVALID_ARG;
        }
    }

    /* Copy A into LU and reject non-finite inputs immediately. */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double v = ML_TENSOR_AT(A, i, j);

            if (ML_UNLIKELY(ml_isnan(v) || ml_isinf(v))) {
                return ML_ERR_NAN_INPUT;
            }

            ML_TENSOR_AT(LU, i, j) = v;
        }

        P[i] = i;
    }

    /* Calculate infinity norm for relative singularity threshold.
     * Scaled accumulation: naive row_sum+=|a| overflows for 1e308
     * entries (valid finite matrix -> spurious NAN_INPUT). Scale by
     * row max first, then rescale. */
    double matrix_norm = 0.0;

    for (int i = 0; i < n; i++) {
        double rmax = 0.0;
        for (int j = 0; j < n; j++) {
            double v = ml_fabs(ML_TENSOR_AT(A, i, j));
            if (v > rmax) rmax = v;
        }
        if (rmax == 0.0) continue;
        {
            double s = 0.0;
            for (int j = 0; j < n; j++) s += ml_fabs(ML_TENSOR_AT(A, i, j)) / rmax;
            {
                double row_sum = s * rmax;
                if (ml_isinf(row_sum)) row_sum = ml_make_inf(0);
                if (row_sum > matrix_norm) matrix_norm = row_sum;
            }
        }
    }

    /*
     * If the norm overflowed, the matrix scale is too large for reliable
     * finite double-precision factorization in this closure path.
     */
    if (ML_UNLIKELY(ml_isinf(matrix_norm))) {
        return ML_ERR_NAN_INPUT;
    }

    /*
     * Exact zero matrix is singular.
     *
     * This must be handled before threshold fallback, otherwise tiny but
     * valid matrices can be falsely classified as singular.
     */
    if (matrix_norm == 0.0) {
        return ML_ERR_SINGULAR;
    }

    /*
     * Relative machine tolerance prevents false singularities on scaled matrices.
     *
     * If the relative threshold underflows to zero, use the smallest positive
     * subnormal instead of an absolute 1e-15 fallback. This preserves scale
     * sensitivity for extremely small matrices.
     */
    double singularity_threshold = matrix_norm * 2.220446049250313e-16 * (double)n;

    if (singularity_threshold == 0.0) {
        singularity_threshold = 4.9406564584124654e-324;
    }

    for (int i = 0; i < n; i++) {
        int max_row = i;
        double max_val = ml_fabs(ML_TENSOR_AT(LU, i, i));

        for (int k = i + 1; k < n; k++) {
            double val = ml_fabs(ML_TENSOR_AT(LU, k, i));

            if (ML_UNLIKELY(ml_isnan(val))) {
                return ML_ERR_SINGULAR;
            }

            if (val > max_val) {
                max_val = val;
                max_row = k;
            }
        }

        if (ML_UNLIKELY(ml_isnan(max_val))) {
            return ML_ERR_SINGULAR;
        }

        if (ML_UNLIKELY(ml_isinf(max_val))) {
            return ML_ERR_NAN_INPUT;
        }

        if (ML_UNLIKELY(max_val < singularity_threshold)) {
            return ML_ERR_SINGULAR;
        }

        if (max_row != i) {
            for (int k = 0; k < n; k++) {
                double tmp = ML_TENSOR_AT(LU, i, k);
                ML_TENSOR_AT(LU, i, k) = ML_TENSOR_AT(LU, max_row, k);
                ML_TENSOR_AT(LU, max_row, k) = tmp;
            }

            int tmp = P[i];
            P[i] = P[max_row];
            P[max_row] = tmp;
        }

        double pivot = ML_TENSOR_AT(LU, i, i);

        if (ML_UNLIKELY(pivot == 0.0 || ml_isnan(pivot))) {
            return ML_ERR_SINGULAR;
        }

        if (ML_UNLIKELY(ml_isinf(pivot))) {
            return ML_ERR_NAN_INPUT;
        }

        for (int k = i + 1; k < n; k++) {
            double mult = ML_TENSOR_AT(LU, k, i) / pivot;
            ML_TENSOR_AT(LU, k, i) = mult;

            for (int j = i + 1; j < n; j++) {
                double updated = ML_TENSOR_AT(LU, k, j) - mult * ML_TENSOR_AT(LU, i, j);

                if (ML_UNLIKELY(!ml_isfinite(updated))) {
                    /* NOTE: ml_types.h defines no ML_ERR_OVERFLOW, so an
                     * overflow-Inf update (|a| > ~1e308 growth) is masked
                     * here as ML_ERR_SINGULAR. Genuine singularity (exact
                     * zero pivot path) and overflow are indistinguishable
                     * at this point without an overflow code. */
                    return ML_ERR_SINGULAR;
                }

                ML_TENSOR_AT(LU, k, j) = updated;
            }
        }
    }

    return ML_SUCCESS;
}

ML_API ml_status_t ml_solve(ml_tensor_view_t A, double* b, double* x, ml_workspace_t* ws) {
    int n = A.rows;

    /* SAFETY BY DEFAULT: Unconditional NULL checks */
    if (ML_UNLIKELY(!A.data || !b || !x || !ws)) {
        return ML_ERR_INVALID_ARG;
    }

    if (ML_UNLIKELY(n <= 0 || A.cols != n)) {
        return ML_ERR_INVALID_ARG;
    }

    /* Reject non-finite RHS entries immediately. */
    for (int i = 0; i < n; i++) {
        if (ML_UNLIKELY(ml_isnan(b[i]) || ml_isinf(b[i]))) {
            return ML_ERR_NAN_INPUT;
        }
    }

    size_t sn = (size_t)n;
    const size_t size_max = (size_t)-1;

    /*
     * Guard against size_t overflow before requesting workspace memory.
     * This is defensive; normal n is far below these limits.
     */
    if (ML_UNLIKELY(sn > size_max / sizeof(double) / sn)) {
        return ML_ERR_WORKSPACE;
    }

    if (ML_UNLIKELY(sn > size_max / sizeof(int))) {
        return ML_ERR_WORKSPACE;
    }

    size_t lu_bytes = sn * sn * sizeof(double);
    size_t p_bytes  = sn * sizeof(int);
    size_t y_bytes  = sn * sizeof(double);

    double* lu_data = (double*)ml_workspace_alloc(ws, lu_bytes);
    int*    P       = (int*)ml_workspace_alloc(ws, p_bytes);
    double* y       = (double*)ml_workspace_alloc(ws, y_bytes);

    if (ML_UNLIKELY(!lu_data || !P || !y)) {
        return ML_ERR_WORKSPACE;
    }

    ml_tensor_view_t LU = ml_tensor_view(lu_data, n, n);

    ml_status_t status = ml_lu_decomp(A, LU, P, ws);
    if (status != ML_SUCCESS) {
        return status;
    }

    /* Forward substitution: Ly = Pb */
    for (int i = 0; i < n; i++) {
        double sum = 0.0;

        for (int j = 0; j < i; j++) {
            sum += ML_TENSOR_AT(LU, i, j) * y[j];
        }

        y[i] = b[P[i]] - sum;

        if (ML_UNLIKELY(!ml_isfinite(y[i]))) {
            return ML_ERR_SINGULAR;
        }
    }

    /* Backward substitution: Ux = y */
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;

        for (int j = i + 1; j < n; j++) {
            sum += ML_TENSOR_AT(LU, i, j) * x[j];
        }

        double pivot = ML_TENSOR_AT(LU, i, i);

        if (ML_UNLIKELY(pivot == 0.0 || ml_isnan(pivot))) {
            return ML_ERR_SINGULAR;
        }

        if (ML_UNLIKELY(ml_isinf(pivot))) {
            return ML_ERR_NAN_INPUT;
        }

        double numerator = y[i] - sum;

        if (ML_UNLIKELY(!ml_isfinite(numerator))) {
            return ML_ERR_SINGULAR;
        }

        x[i] = numerator / pivot;

        if (ML_UNLIKELY(!ml_isfinite(x[i]))) {
            return ML_ERR_SINGULAR;
        }
    }

    return ML_SUCCESS;
}

/* Matrix-Vector Multiplication (y = Ax) */
ML_API void ml_matvec(ml_tensor_view_t A, const double* x, double* out) {
    int n = A.rows;
    int m = A.cols;

    if (ML_UNLIKELY(!A.data || !x || !out)) {
        return;
    }

    if (ML_UNLIKELY(n <= 0 || m <= 0)) {
        return;
    }

    for (int i = 0; i < n; i++) {
        double sum = 0.0;

        for (int j = 0; j < m; j++) {
            sum += ML_TENSOR_AT(A, i, j) * x[j];
        }

        out[i] = sum;
    }
}

ML_API ml_status_t ml_cholesky(ml_tensor_view_t A, ml_tensor_view_t L) {
    int n = A.rows;
    if (ML_UNLIKELY(n <= 0 || !A.data || !L.data)) return ML_ERR_INVALID_ARG;
    if (ML_UNLIKELY(A.cols != n || L.rows != n || L.cols != n)) return ML_ERR_INVALID_ARG;
    /* Symmetry check: Cholesky requires symmetric input; previously the
     * upper triangle was silently ignored. Reject asymmetric matrices. */
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double u = ML_TENSOR_AT(A, i, j);
            double l = ML_TENSOR_AT(A, j, i);
            if (ML_UNLIKELY(!ml_isfinite(u) || !ml_isfinite(l))) return ML_ERR_NAN_INPUT;
            {
                double diff = ml_fabs(u - l);
                double sc = ml_fabs(u) > ml_fabs(l) ? ml_fabs(u) : ml_fabs(l);
                if (sc == 0.0) continue;
                if (!(diff <= 1e-12 * sc)) return ML_ERR_INVALID_ARG;
            }
        }
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            double s = ML_TENSOR_AT(A, i, j);
            if (ML_UNLIKELY(!ml_isfinite(s))) return ML_ERR_NAN_INPUT;
            for (int k = 0; k < j; k++) {
                s -= ML_TENSOR_AT(L, i, k) * ML_TENSOR_AT(L, j, k);
            }
            if (i == j) {
                if (ML_UNLIKELY(!(s > 0.0) || !ml_isfinite(s))) return ML_ERR_SINGULAR;
                ML_TENSOR_AT(L, i, j) = ml_sqrt(s);
            } else {
                double d = ML_TENSOR_AT(L, j, j);
                if (ML_UNLIKELY(d == 0.0 || !ml_isfinite(d))) return ML_ERR_SINGULAR;
                ML_TENSOR_AT(L, i, j) = s / d;
            }
        }
        for (int j = i + 1; j < n; j++) {
            ML_TENSOR_AT(L, i, j) = 0.0;
        }
    }
    return ML_SUCCESS;
}

ML_API ml_status_t ml_qr_solve(ml_tensor_view_t A, const double* b, double* x, ml_workspace_t* ws) {
    /* Modified Gram-Schmidt least-squares: works for square nonsingular
     * and tall m>=n full-rank. Q stored implicitly in V, R upper. */
    int m = A.rows, n = A.cols;
    if (ML_UNLIKELY(m <= 0 || n <= 0 || m < n)) return ML_ERR_INVALID_ARG;
    if (ML_UNLIKELY(!A.data || !b || !x || !ws)) return ML_ERR_INVALID_ARG;
    for (int i = 0; i < m; i++) {
        if (ML_UNLIKELY(!ml_isfinite(b[i]))) return ML_ERR_NAN_INPUT;
    }
    {
        size_t sm = (size_t)m, sn = (size_t)n;
        const size_t size_max = (size_t)-1;
        /* size_t overflow guards (same style as ml_solve): reject
         * sm*sn / sn*sn element counts that would wrap before alloc. */
        if (ML_UNLIKELY(sn > 0 && sm > size_max / sizeof(double) / sn)) {
            return ML_ERR_WORKSPACE;
        }
        if (ML_UNLIKELY(sn > size_max / sizeof(double) / sn)) {
            return ML_ERR_WORKSPACE;
        }
        if (ML_UNLIKELY(sn > size_max / sizeof(double))) {
            return ML_ERR_WORKSPACE;
        }
        size_t v_bytes = sm * sn * sizeof(double);
        size_t r_bytes = sn * sn * sizeof(double);
        size_t q = sn * sizeof(double);
        double* V = (double*)ml_workspace_alloc(ws, v_bytes);
        double* R = (double*)ml_workspace_alloc(ws, r_bytes);
        double* y = (double*)ml_workspace_alloc(ws, q);
        if (ML_UNLIKELY(!V || !R || !y)) return ML_ERR_WORKSPACE;
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                double v = ML_TENSOR_AT(A, i, j);
                if (ML_UNLIKELY(!ml_isfinite(v))) return ML_ERR_NAN_INPUT;
                V[(size_t)i * sn + (size_t)j] = v;
            }
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) R[(size_t)i * sn + (size_t)j] = 0.0;
        }
        for (int k = 0; k < n; k++) {
            /* Max-scaled column norm: naive sum v*v overflows for 1e200
             * columns though the QR factor exists. */
            double cmax = 0.0;
            for (int i = 0; i < m; i++) {
                double av = ml_fabs(V[(size_t)i * sn + (size_t)k]);
                if (av > cmax) cmax = av;
            }
            if (!(cmax > 0.0) || !ml_isfinite(cmax)) return ML_ERR_SINGULAR;
            {
                double s = 0.0;
                for (int i = 0; i < m; i++) {
                    double r = V[(size_t)i * sn + (size_t)k] / cmax;
                    s += r * r;
                }
                {
                    double nrm = cmax * ml_sqrt(s);
                    if (ML_UNLIKELY(!(nrm > 0.0) || !ml_isfinite(nrm))) return ML_ERR_SINGULAR;
                    R[(size_t)k * sn + (size_t)k] = nrm;
                    for (int i = 0; i < m; i++) V[(size_t)i * sn + (size_t)k] /= nrm;
                }
            }
            for (int j = k + 1; j < n; j++) {
                double dot = 0.0;
                for (int i = 0; i < m; i++) dot += V[(size_t)i * sn + (size_t)k] * V[(size_t)i * sn + (size_t)j];
                R[(size_t)k * sn + (size_t)j] = dot;
                for (int i = 0; i < m; i++) V[(size_t)i * sn + (size_t)j] -= dot * V[(size_t)i * sn + (size_t)k];
            }
        }
        for (int k = 0; k < n; k++) {
            double dot = 0.0;
            for (int i = 0; i < m; i++) dot += V[(size_t)i * sn + (size_t)k] * b[i];
            y[k] = dot;
        }
        for (int i = n - 1; i >= 0; i--) {
            double s = y[i];
            for (int j = i + 1; j < n; j++) s -= R[(size_t)i * sn + (size_t)j] * x[j];
            {
                double d = R[(size_t)i * sn + (size_t)i];
                if (ML_UNLIKELY(d == 0.0 || !ml_isfinite(d))) return ML_ERR_SINGULAR;
                x[i] = s / d;
                if (ML_UNLIKELY(!ml_isfinite(x[i]))) return ML_ERR_SINGULAR;
            }
        }
        return ML_SUCCESS;
    }
}

ML_API ml_status_t ml_solve_refined(ml_tensor_view_t A, double* b, double* x, ml_workspace_t* ws) {
    /* One step of iterative refinement: x0=solve, r=b-Ax0 (Kahan),
     * solve A*d=r, x=x0+d. Halves forward error for ill-conditioned.
     * Workspace: ~2x one ml_solve (two extra n-vectors plus inner solve). */
    ml_status_t st = ml_solve(A, b, x, ws);
    if (st != ML_SUCCESS) return st;
    {
        int n = A.rows;
        size_t sn = (size_t)n;
        const size_t size_max = (size_t)-1;
        /* size_t overflow guard (same style as ml_solve). */
        if (ML_UNLIKELY(sn > size_max / sizeof(double))) {
            return ML_ERR_WORKSPACE;
        }
        double* r = (double*)ml_workspace_alloc(ws, sn * sizeof(double));
        double* d = (double*)ml_workspace_alloc(ws, sn * sizeof(double));
        if (ML_UNLIKELY(!r || !d)) return ML_ERR_WORKSPACE;
        for (int i = 0; i < n; i++) {
            double s = b[i];
            double comp = 0.0;
            for (int j = 0; j < n; j++) {
                double p = -ML_TENSOR_AT(A, i, j) * x[j];
                double yy = p - comp;
                double t = s + yy;
                comp = (t - s) - yy;
                s = t;
            }
            r[i] = s;
            if (!ml_isfinite(r[i])) return st;
        }
        {
            double rn = 0.0, bn = 0.0;
            for (int i = 0; i < n; i++) {
                rn += r[i] * r[i];
                bn += b[i] * b[i];
            }
            if (!(rn < 1e-28 * bn + 1e-300)) {
                ml_tensor_view_t A2 = A;
                (void)A2;
                /* Solve A*d=r via fresh LU (uses its own workspace allocs). */
                /* Build temporary views: need square solve, so call ml_solve
                 * on a copy? ml_solve takes A + rhs; reuse A directly. */
                ml_status_t st2 = ml_solve(A, r, d, ws);
                if (st2 != ML_SUCCESS) return st;
                for (int i = 0; i < n; i++) {
                    double nx = x[i] + d[i];
                    if (ml_isfinite(nx)) x[i] = nx;
                }
            }
        }
        return ML_SUCCESS;
    }
}

ML_API double ml_determinant(ml_tensor_view_t A, ml_workspace_t* ws) {
    int n = A.rows;
    if (ML_UNLIKELY(n <= 0 || !A.data || !ws)) return ml_make_nan();
    if (ML_UNLIKELY(A.cols != n)) return ml_make_nan();
    {
        size_t sn = (size_t)n;
        const size_t size_max = (size_t)-1;
        /* size_t overflow guards (same style as ml_solve). */
        if (ML_UNLIKELY(sn > size_max / sizeof(double) / sn)) return ml_make_nan();
        if (ML_UNLIKELY(sn > size_max / sizeof(int))) return ml_make_nan();
        double* lu = (double*)ml_workspace_alloc(ws, sn * sn * sizeof(double));
        int* P = (int*)ml_workspace_alloc(ws, sn * sizeof(int));
        if (ML_UNLIKELY(!lu || !P)) return ml_make_nan();
        {
            ml_tensor_view_t LU = ml_tensor_view(lu, n, n);
            ml_status_t st = ml_lu_decomp(A, LU, P, ws);
            if (st != ML_SUCCESS) return (st == ML_ERR_SINGULAR) ? 0.0 : ml_make_nan();
            {
                double det = 1.0;
                int swaps = 0;
                for (int i = 0; i < n; i++) {
                    if (P[i] != i) swaps++;
                    det *= ML_TENSOR_AT(LU, i, i);
                }
                /* Parity via permutation sign: count inversions of P. */
                {
                    int inv = 0;
                    for (int i = 0; i < n; i++) {
                        for (int j = i + 1; j < n; j++) {
                            if (P[i] > P[j]) inv++;
                        }
                    }
                    if ((inv % 2) != 0) det = -det;
                    (void)swaps;
                }
                return det;
            }
        }
    }
}

ML_API ml_status_t ml_inverse(ml_tensor_view_t A, ml_tensor_view_t Inv, ml_workspace_t* ws) {
    int n = A.rows;
    if (ML_UNLIKELY(n <= 0 || !A.data || !Inv.data || !ws)) return ML_ERR_INVALID_ARG;
    if (ML_UNLIKELY(A.cols != n || Inv.rows != n || Inv.cols != n)) return ML_ERR_INVALID_ARG;
    {
        size_t sn = (size_t)n;
        const size_t size_max = (size_t)-1;
        /* size_t overflow guards (same style as ml_solve). */
        if (ML_UNLIKELY(sn > size_max / sizeof(double) / sn)) {
            return ML_ERR_WORKSPACE;
        }
        if (ML_UNLIKELY(sn > size_max / sizeof(int))) {
            return ML_ERR_WORKSPACE;
        }
        if (ML_UNLIKELY(sn > size_max / sizeof(double))) {
            return ML_ERR_WORKSPACE;
        }
        double* lu = (double*)ml_workspace_alloc(ws, sn * sn * sizeof(double));
        int* P = (int*)ml_workspace_alloc(ws, sn * sizeof(int));
        double* e = (double*)ml_workspace_alloc(ws, sn * sizeof(double));
        double* col = (double*)ml_workspace_alloc(ws, sn * sizeof(double));
        if (ML_UNLIKELY(!lu || !P || !e || !col)) return ML_ERR_WORKSPACE;
        {
            ml_tensor_view_t LU = ml_tensor_view(lu, n, n);
            ml_status_t st = ml_lu_decomp(A, LU, P, ws);
            if (st != ML_SUCCESS) return st;
            for (int j = 0; j < n; j++) {
                for (int i = 0; i < n; i++) e[i] = (P[i] == j) ? 1.0 : 0.0;
                /* Forward: Ly=Pb(=e since e already permuted). */
                for (int i = 0; i < n; i++) {
                    double s = e[i];
                    for (int k = 0; k < i; k++) s -= ML_TENSOR_AT(LU, i, k) * col[k];
                    col[i] = s;
                    if (ML_UNLIKELY(!ml_isfinite(col[i]))) return ML_ERR_SINGULAR;
                }
                for (int i = n - 1; i >= 0; i--) {
                    double s = col[i];
                    for (int k = i + 1; k < n; k++) s -= ML_TENSOR_AT(LU, i, k) * ML_TENSOR_AT(Inv, k, j);
                    {
                        double d = ML_TENSOR_AT(LU, i, i);
                        if (ML_UNLIKELY(d == 0.0 || !ml_isfinite(d))) return ML_ERR_SINGULAR;
                        ML_TENSOR_AT(Inv, i, j) = s / d;
                    }
                }
            }
            return ML_SUCCESS;
        }
    }
}

ML_API void ml_transpose(ml_tensor_view_t A, ml_tensor_view_t T) {
    if (ML_UNLIKELY(!A.data || !T.data)) return;
    if (ML_UNLIKELY(T.rows != A.cols || T.cols != A.rows)) return;
    {
        int m = A.rows, n = A.cols;
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                ML_TENSOR_AT(T, j, i) = ML_TENSOR_AT(A, i, j);
            }
        }
    }
}

ML_API double ml_dot(const double* a, const double* b, int n) {
    if (ML_UNLIKELY(!a || !b || n <= 0)) return ml_make_nan();
    {
        double s = 0.0, c = 0.0;
        for (int i = 0; i < n; i++) {
            if (ML_UNLIKELY(!ml_isfinite(a[i]) || !ml_isfinite(b[i]))) return ml_make_nan();
            double p = a[i] * b[i] - c;
            double t = s + p;
            c = (t - s) - p;
            s = t;
        }
        return s;
    }
}

ML_API void ml_cross3(const double a[3], const double b[3], double out[3]) {
    if (ML_UNLIKELY(!a || !b || !out)) return;
    out[0] = a[1] * b[2] - a[2] * b[1];
    out[1] = a[2] * b[0] - a[0] * b[2];
    out[2] = a[0] * b[1] - a[1] * b[0];
}

ML_API double ml_norm(const double* a, int n) {
    if (ML_UNLIKELY(!a || n <= 0)) return ml_make_nan();
    {
        double m = 0.0;
        for (int i = 0; i < n; i++) {
            if (ML_UNLIKELY(!ml_isfinite(a[i]))) return ml_make_nan();
            double v = ml_fabs(a[i]);
            if (v > m) m = v;
        }
        if (m == 0.0) return 0.0;
        {
            double s = 0.0;
            for (int i = 0; i < n; i++) {
                double r = a[i] / m;
                s += r * r;
            }
            return m * ml_sqrt(s);
        }
    }
}

ML_API double ml_point_line_dist2d(double px, double py, double ax, double ay, double bx, double by) {
    if (ML_UNLIKELY(!ml_isfinite(px) || !ml_isfinite(py) || !ml_isfinite(ax) || !ml_isfinite(ay) || !ml_isfinite(bx) || !ml_isfinite(by))) {
        return ml_make_nan();
    }
    {
        double dx = bx - ax, dy = by - ay;
        double L = ml_hypot_internal(dx, dy);
        if (L == 0.0 || !ml_isfinite(L)) return ml_make_nan();
        {
            /* Numerator dx*(ay-py)-dy*(ax-px) overflows for 1e200 coords
             * though the ratio is finite; evaluate in long double. */
            long double num = (long double)dx * ((long double)ay - (long double)py)
                            - (long double)dy * ((long double)ax - (long double)px);
            if (num < 0.0L) num = -num;
            long double r = num / (long double)L;
            if (r > (long double)1.7976931348623157e308L) return ml_make_inf(0);
            return (double)r;
        }
    }
}

ML_API double ml_point_plane_dist(double px, double py, double pz, double a, double b, double c, double d) {
    if (ML_UNLIKELY(!ml_isfinite(px) || !ml_isfinite(py) || !ml_isfinite(pz) || !ml_isfinite(a) || !ml_isfinite(b) || !ml_isfinite(c) || !ml_isfinite(d))) {
        return ml_make_nan();
    }
    {
        /* Scaled 3-norm: naive sqrt(a^2+b^2+c^2) overflows for 1e200. */
        double m = ml_fabs(a);
        if (ml_fabs(b) > m) m = ml_fabs(b);
        if (ml_fabs(c) > m) m = ml_fabs(c);
        if (!(m > 0.0) || !ml_isfinite(m)) return ml_make_nan();
        {
            double ra = a / m, rb = b / m, rc = c / m;
            double n = m * ml_sqrt(ra*ra + rb*rb + rc*rc);
            if (n == 0.0 || !ml_isfinite(n)) return ml_make_nan();
            {
                long double num = (long double)a * (long double)px + (long double)b * (long double)py
                                + (long double)c * (long double)pz + (long double)d;
                if (num < 0.0L) num = -num;
                long double r = num / (long double)n;
                if (r > (long double)1.7976931348623157e308L) return ml_make_inf(0);
                return (double)r;
            }
        }
    }
}

ML_API ml_status_t ml_eigen2x2(double a, double b, double c, double d, double *l0, double *l1) {
    if (ML_UNLIKELY(!l0 || !l1)) return ML_ERR_INVALID_ARG;
    if (ML_UNLIKELY(!ml_isfinite(a) || !ml_isfinite(b) || !ml_isfinite(c) || !ml_isfinite(d))) {
        *l0 = ml_make_nan(); *l1 = ml_make_nan(); return ML_ERR_NAN_INPUT;
    }
    {
        /* Scaled trace/det: naive tr*tr-4det overflows for 1e200 diag
         * and (tr-s)/2 cancels when |tr|>>sqrt|det|. Use long double
         * discriminant plus stable q-form. */
        long double tr = (long double)a + (long double)d;
        long double det = (long double)a * (long double)d - (long double)b * (long double)c;
        long double disc = tr * tr - 4.0L * det;
        if (!(disc >= 0.0L) || !(disc < (long double)1e4932L)) {
            *l0 = ml_make_nan(); *l1 = ml_make_nan(); return ML_ERR_SINGULAR;
        }
        {
            double s = (double)__builtin_sqrtl(disc);
            if (!ml_isfinite(s)) { *l0 = ml_make_nan(); *l1 = ml_make_nan(); return ML_ERR_SINGULAR; }
            {
                double trd = a + d;
                if (!ml_isfinite(trd)) { *l0 = ml_make_nan(); *l1 = ml_make_nan(); return ML_ERR_SINGULAR; }
                /* Stable: q=0.5*(tr+sgn(tr)*s) is the large-magnitude root
                 * (no cancellation); other root = det/q in long double
                 * so det=1e400 (overflowing double) still divides
                 * correctly to 1e200. */
                double q = 0.5 * (trd + (trd >= 0.0 ? s : -s));
                if (q == 0.0) {
                    *l0 = (trd + s) * 0.5;
                    *l1 = (trd - s) * 0.5;
                } else {
                    long double l1l = det / (long double)q;
                    *l0 = q;
                    if (l1l > (long double)1.7976931348623157e308L) *l1 = ml_make_inf(0);
                    else if (l1l < -(long double)1.7976931348623157e308L) *l1 = ml_make_inf(1);
                    else *l1 = (double)l1l;
                    /* det/q would overflow if q tiny; fall back to sum form. */
                    if (!ml_isfinite(*l1)) *l1 = (trd - s) * 0.5;
                    /* Order l0>=l1 for determinism when tr>=0, else keep
                     * q-form pair (product still = det). */
                    if (*l1 > *l0) { double t = *l0; *l0 = *l1; *l1 = t; }
                }
                return ML_SUCCESS;
            }
        }
    }
}

ML_API ml_status_t ml_jacobi_eigen_symmetric(ml_tensor_view_t A, double *evals,
                                             ml_tensor_view_t V, int max_sweeps) {
    /* Cyclic Jacobi: spectral theorem A=Q diag Q^T. Backward stable;
     * off(A) decreases monotonically; ~5-10 sweeps for n<=32. */
    int n = A.rows;
    if (ML_UNLIKELY(n <= 0 || !A.data || !evals)) return ML_ERR_INVALID_ARG;
    if (ML_UNLIKELY(A.cols != n)) return ML_ERR_INVALID_ARG;
    if (ML_UNLIKELY(max_sweeps <= 0)) max_sweeps = 20;
    if (ML_UNLIKELY(max_sweeps > 100)) max_sweeps = 100;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double v = ML_TENSOR_AT(A, i, j);
            if (ML_UNLIKELY(!ml_isfinite(v))) return ML_ERR_NAN_INPUT;
        }
        if (V.data) {
            for (int j = 0; j < n; j++) {
                if (V.rows != n || V.cols != n) return ML_ERR_INVALID_ARG;
                ML_TENSOR_AT(V, i, j) = (i == j) ? 1.0 : 0.0;
            }
        }
    }
    /* Symmetry check (relative 1e-12). */
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double u = ML_TENSOR_AT(A, i, j), l = ML_TENSOR_AT(A, j, i);
            double sc = ml_fabs(u) > ml_fabs(l) ? ml_fabs(u) : ml_fabs(l);
            if (sc == 0.0) continue;
            if (!(ml_fabs(u - l) <= 1e-12 * sc)) return ML_ERR_INVALID_ARG;
        }
    }
    {
        size_t sn = (size_t)n;
        /* Stack-local work copy (thread-safe; no static storage).
         * Capacity 4096 doubles => n<=64 cap (64*64 == 4096). */
        double W[1024 * 4];
        double *w = W;
        int need = n * n;
        if (need > (int)(sizeof(W) / sizeof(W[0]))) return ML_ERR_WORKSPACE;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) w[(size_t)i * sn + (size_t)j] = ML_TENSOR_AT(A, i, j);
        }
        /* Scale-invariant stop: off <= eps^2 * ||A||_F^2, with the squared
         * Frobenius norm captured at entry (long double accumulation). */
        long double fnorm = 0.0L;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                long double v = (long double)w[(size_t)i * sn + (size_t)j];
                fnorm += v * v;
            }
        }
        {
            const long double eps = 2.220446049250313e-16L;
            const long double tol = eps * eps * fnorm;
            for (int sw = 0; sw < max_sweeps; sw++) {
            long double off = 0.0L;
            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    long double v = (long double)w[(size_t)i * sn + (size_t)j];
                    off += v * v;
                }
            }
            if (!(off > 0.0L)) break;
            if (off <= tol) break;
            for (int p = 0; p < n - 1; p++) {
                for (int q2 = p + 1; q2 < n; q2++) {
                    double apq = w[(size_t)p * sn + (size_t)q2];
                    if (apq == 0.0) continue;
                    double app = w[(size_t)p * sn + (size_t)p];
                    double aqq = w[(size_t)q2 * sn + (size_t)q2];
                    double theta = (aqq - app) / (2.0 * apq);
                    double t = ((theta >= 0.0) ? 1.0 : -1.0)
                        / (ml_fabs(theta) + ml_sqrt(theta * theta + 1.0));
                    double c = 1.0 / ml_sqrt(t * t + 1.0);
                    double s = t * c;
                    for (int k = 0; k < n; k++) {
                        if (k != p && k != q2) {
                            double akp = w[(size_t)k * sn + (size_t)p];
                            double akq = w[(size_t)k * sn + (size_t)q2];
                            w[(size_t)k * sn + (size_t)p] = c * akp - s * akq;
                            w[(size_t)p * sn + (size_t)k] = w[(size_t)k * sn + (size_t)p];
                            w[(size_t)k * sn + (size_t)q2] = s * akp + c * akq;
                            w[(size_t)q2 * sn + (size_t)k] = w[(size_t)k * sn + (size_t)q2];
                        }
                    }
                    {
                        double npp = c * c * app - 2.0 * s * c * apq + s * s * aqq;
                        double nqq = s * s * app + 2.0 * s * c * apq + c * c * aqq;
                        w[(size_t)p * sn + (size_t)p] = npp;
                        w[(size_t)q2 * sn + (size_t)q2] = nqq;
                        w[(size_t)p * sn + (size_t)q2] = 0.0;
                        w[(size_t)q2 * sn + (size_t)p] = 0.0;
                    }
                    if (V.data) {
                        for (int k = 0; k < n; k++) {
                            double vkp = ML_TENSOR_AT(V, k, p);
                            double vkq = ML_TENSOR_AT(V, k, q2);
                            ML_TENSOR_AT(V, k, p) = c * vkp - s * vkq;
                            ML_TENSOR_AT(V, k, q2) = s * vkp + c * vkq;
                        }
                    }
                }
            }
            }
        }
        for (int i = 0; i < n; i++) evals[i] = w[(size_t)i * sn + (size_t)i];
        return ML_SUCCESS;
    }
}

ML_API ml_status_t ml_svd_2x2(double a, double b, double c, double d,
                              double *s0, double *s1) {
    if (ML_UNLIKELY(!s0 || !s1)) return ML_ERR_INVALID_ARG;
    if (ML_UNLIKELY(!ml_isfinite(a) || !ml_isfinite(b) || !ml_isfinite(c) || !ml_isfinite(d))) {
        *s0 = ml_make_nan(); *s1 = ml_make_nan(); return ML_ERR_NAN_INPUT;
    }
    {
        /* s = sqrt(eig(A^T A)); A^T A symmetric PSD: use stable eigen. */
        long double e = (long double)a * (long double)a + (long double)c * (long double)c;
        long double f = (long double)a * (long double)b + (long double)c * (long double)d;
        long double g = (long double)b * (long double)b + (long double)d * (long double)d;
        long double tr = e + g;
        long double det = e * g - f * f;
        if (det < 0.0L && det > -1e-18L * (e * g + 1.0L)) det = 0.0L;
        if (!(det >= 0.0L)) { *s0 = ml_make_nan(); *s1 = ml_make_nan(); return ML_ERR_SINGULAR; }
        {
            long double disc = tr * tr - 4.0L * det;
            if (disc < 0.0L && disc > -1e-18L * tr * tr) disc = 0.0L;
            if (!(disc >= 0.0L)) { *s0 = ml_make_nan(); *s1 = ml_make_nan(); return ML_ERR_SINGULAR; }
            long double sq = __builtin_sqrtl(disc);
            long double l0 = (tr + sq) * 0.5L, l1 = (tr - sq) * 0.5L;
            if (l1 < 0.0L && l1 > -1e-18L * tr) l1 = 0.0L;
            if (l0 < 0.0L || l1 < 0.0L) { *s0 = ml_make_nan(); *s1 = ml_make_nan(); return ML_ERR_SINGULAR; }
            *s0 = (double)__builtin_sqrtl(l0);
            *s1 = (double)__builtin_sqrtl(l1);
            if (!ml_isfinite(*s0) || !ml_isfinite(*s1)) return ML_ERR_SINGULAR;
            if (*s1 > *s0) { double t = *s0; *s0 = *s1; *s1 = t; }
            return ML_SUCCESS;
        }
    }
}

ML_API ml_status_t ml_matrix_exp_2x2(double a, double b, double c, double d,
                                     double *e00, double *e01,
                                     double *e10, double *e11) {
    if (ML_UNLIKELY(!e00 || !e01 || !e10 || !e11)) return ML_ERR_INVALID_ARG;
    if (ML_UNLIKELY(!ml_isfinite(a) || !ml_isfinite(b) || !ml_isfinite(c) || !ml_isfinite(d))) return ML_ERR_NAN_INPUT;
    {
        /* Closed form via trace/discriminant (Jordan/cayley-hamilton):
         * m=(a+d)/2, D=((a-d)/2)^2+bc, exp(A)=e^m(cosh+...). */
        long double m = ((long double)a + (long double)d) * 0.5L;
        long double ah = ((long double)a - (long double)d) * 0.5L;
        long double D = ah * ah + (long double)b * (long double)c;
        /* e^m computed separately so the overflow path can inspect it
         * per element below. */
        long double em = __builtin_expl(m);
        long double l00 = 0.0L, l01 = 0.0L, l10 = 0.0L, l11 = 0.0L;
        if (D >= 0.0L) {
            long double s = __builtin_sqrtl(D);
            long double ch = (s == 0.0L) ? 1.0L : __builtin_coshl(s);
            long double sh = (s == 0.0L) ? 1.0L : __builtin_sinhl(s) / s;
            /* Avoid 0*Inf -> NaN when ah==0 and sh==Inf (or b/c==0):
             * the exact term is 0 in those cases. */
            long double ah_sh = (ah == 0.0L) ? 0.0L : ah * sh;
            l00 = em * (ch + ah_sh);
            l01 = (b == 0.0) ? 0.0L : em * (long double)b * sh;
            l10 = (c == 0.0) ? 0.0L : em * (long double)c * sh;
            l11 = em * (ch - ah_sh);
            *e00 = (double)l00;
            *e01 = (double)l01;
            *e10 = (double)l10;
            *e11 = (double)l11;
        } else {
            long double s = __builtin_sqrtl(-D);
            long double cs = __builtin_cosl(s);
            long double sn = (s == 0.0L) ? 1.0L : __builtin_sinl(s) / s;
            long double ah_sn = (ah == 0.0L) ? 0.0L : ah * sn;
            l00 = em * (cs + ah_sn);
            l01 = (b == 0.0) ? 0.0L : em * (long double)b * sn;
            l10 = (c == 0.0) ? 0.0L : em * (long double)c * sn;
            l11 = em * (cs - ah_sn);
            *e00 = (double)l00;
            *e01 = (double)l01;
            *e10 = (double)l10;
            *e11 = (double)l11;
        }
        if (!ml_isfinite(*e00) || !ml_isfinite(*e01) || !ml_isfinite(*e10) || !ml_isfinite(*e11)) {
            /* Per-element overflow fixup: diagonal entries round to signed
             * Inf; off-diagonal em*b*sh with b==0 (resp. c==0) is exactly
             * 0, not Inf/NaN, even when em is +Inf. */
            if (!ml_isfinite(*e00)) {
                if (l00 == 0.0L) {
                    *e00 = 0.0;
                } else if (l00 != l00) {
                    return ML_ERR_SINGULAR;
                } else {
                    *e00 = ml_make_inf((l00 < 0.0L) ? 1 : 0);
                }
            }
            if (!ml_isfinite(*e01)) {
                if (b == 0.0) {
                    *e01 = 0.0;
                } else if (l01 != l01) {
                    return ML_ERR_SINGULAR;
                } else if (l01 == 0.0L) {
                    *e01 = 0.0;
                } else {
                    *e01 = ml_make_inf((l01 < 0.0L) ? 1 : 0);
                }
            }
            if (!ml_isfinite(*e10)) {
                if (c == 0.0) {
                    *e10 = 0.0;
                } else if (l10 != l10) {
                    return ML_ERR_SINGULAR;
                } else if (l10 == 0.0L) {
                    *e10 = 0.0;
                } else {
                    *e10 = ml_make_inf((l10 < 0.0L) ? 1 : 0);
                }
            }
            if (!ml_isfinite(*e11)) {
                if (l11 == 0.0L) {
                    *e11 = 0.0;
                } else if (l11 != l11) {
                    return ML_ERR_SINGULAR;
                } else {
                    *e11 = ml_make_inf((l11 < 0.0L) ? 1 : 0);
                }
            }
            /* Overflow rounded to signed Inf is representable: only NaN
             * remains a failure here. */
            if (ml_isnan(*e00) || ml_isnan(*e01) ||
                ml_isnan(*e10) || ml_isnan(*e11)) {
                return ML_ERR_SINGULAR;
            }
            return ML_SUCCESS;
        }
        return ML_SUCCESS;
    }
}
