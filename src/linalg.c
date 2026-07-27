#include "ml_compiler.h"
#include "ml_linalg.h"

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

    /* Calculate infinity norm for relative singularity threshold */
    double matrix_norm = 0.0;

    for (int i = 0; i < n; i++) {
        double row_sum = 0.0;

        for (int j = 0; j < n; j++) {
            row_sum += ml_fabs(ML_TENSOR_AT(A, i, j));
        }

        if (row_sum > matrix_norm) {
            matrix_norm = row_sum;
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

                if (ML_UNLIKELY(ml_isnan(updated))) {
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
