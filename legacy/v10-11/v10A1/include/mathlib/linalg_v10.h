#ifndef MATHLIB_V10_LINALG_H
#define MATHLIB_V10_LINALG_H

#include "tensor.h"

// Zero-Allocation LU Decomposition with Partial Pivoting
// L and U are packed into a single workspace buffer to save memory.
// P is the permutation vector.
static inline int ml_lu_decomp_v10(ml_tensor_view_t A, ml_tensor_view_t LU, int* P, ml_workspace_t* ws) {
    (void)ws; // Suppress unused warning
    int n = A.rows;

    // Copy A into LU
    for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
            ML_TENSOR_AT(LU, i, j) = ML_TENSOR_AT(A, i, j);
        }
        P[i] = i;
    }

    for (int i = 0; i < n; i++) {
        // Partial Pivoting
        int max_row = i;
        double max_val = ML_TENSOR_AT(LU, i, i);
        if (max_val < 0) max_val = -max_val;

        for (int k = i + 1; k < n; k++) {
            double val = ML_TENSOR_AT(LU, k, i);
            if (val < 0) val = -val;
            if (val > max_val) {
                max_val = val;
                max_row = k;
            }
        }
        if (max_val == 0.0) return -1; // Singular

        // Swap rows in LU and P
        if (max_row != i) {
            for (int k = 0; k < n; k++) {
                double tmp = ML_TENSOR_AT(LU, i, k);
                ML_TENSOR_AT(LU, i, k) = ML_TENSOR_AT(LU, max_row, k);
                ML_TENSOR_AT(LU, max_row, k) = tmp;
            }
            int tmp = P[i]; P[i] = P[max_row]; P[max_row] = tmp;
        }

        // Elimination (Store multipliers in the lower triangle of LU)
        double pivot = ML_TENSOR_AT(LU, i, i);
        for (int k = i + 1; k < n; k++) {
            double mult = ML_TENSOR_AT(LU, k, i) / pivot;
            ML_TENSOR_AT(LU, k, i) = mult; // Store L
            for (int j = i + 1; j < n; j++) {
                ML_TENSOR_AT(LU, k, j) -= mult * ML_TENSOR_AT(LU, i, j);
            }
        }
    }
    return 0;
}

// Zero-Allocation Linear Solve (Ax = b)
static inline int ml_solve_v10(ml_tensor_view_t A, double* b, double* x, ml_workspace_t* ws) {
    int n = A.rows;

    // Allocate scratchpad from workspace (NO MALLOC)
    double* lu_data = (double*)ml_workspace_alloc(ws, n * n * sizeof(double));
    int* P = (int*)ml_workspace_alloc(ws, n * sizeof(int));
    double* y = (double*)ml_workspace_alloc(ws, n * sizeof(double));

    if (!lu_data || !P || !y) return -2; // Workspace too small

    ml_tensor_view_t LU = ml_tensor_view(lu_data, n, n);

    if (ml_lu_decomp_v10(A, LU, P, ws) != 0) return -1;

    // Forward substitution (Ly = Pb)
    for (int i = 0; i < n; i++) {
        double sum = 0;
        for (int j = 0; j < i; j++) sum += ML_TENSOR_AT(LU, i, j) * y[j];
        y[i] = b[P[i]] - sum;
    }

    // Backward substitution (Ux = y)
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0;
        for (int j = i + 1; j < n; j++) sum += ML_TENSOR_AT(LU, i, j) * x[j];
        x[i] = (y[i] - sum) / ML_TENSOR_AT(LU, i, i);
    }

    return 0;
}
#endif
