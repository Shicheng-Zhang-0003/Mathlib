#ifndef LIBMATHC_MATRIX_H
#define LIBMATHC_MATRIX_H

#include <stdlib.h>
#include <string.h>

typedef struct {
    int rows, cols;
    double *data;
} ml_matrix;

#define ML_MAT(m, r, c) ((m).data[(r) * (m).cols + (c)])

inline ml_matrix ml_matrix_create(int r, int c) {
    ml_matrix m;
    m.rows = r; m.cols = c;
    m.data = (double*)calloc(r * c, sizeof(double));
    return m;
}

inline void ml_matrix_free(ml_matrix m) {
    if(m.data) free(m.data);
}

inline ml_matrix ml_matrix_identity(int n) {
    ml_matrix m = ml_matrix_create(n, n);
    for (int i = 0; i < n; i++) ML_MAT(m, i, i) = 1.0;
    return m;
}

inline ml_matrix ml_matrix_mul(ml_matrix a, ml_matrix b) {
    ml_matrix res = ml_matrix_create(a.rows, b.cols);
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < b.cols; j++) {
            double sum = 0;
            for (int k = 0; k < a.cols; k++) {
                sum += ML_MAT(a, i, k) * ML_MAT(b, k, j);
            }
            ML_MAT(res, i, j) = sum;
        }
    }
    return res;
}


// --- Advanced NxN Decompositions & Solvers ---

// LU Decomposition with Partial Pivoting (PA = LU)
// Returns 0 on success, -1 if singular.
inline int ml_matrix_lu(ml_matrix A, ml_matrix *L, ml_matrix *U, int *P) {
    int n = A.rows;
    *L = ml_matrix_create(n, n);
    *U = ml_matrix_create(n, n);
    for(int i=0; i<n; i++) { P[i] = i; ML_MAT(*L, i, i) = 1.0; }

    // Copy A to U
    for(int i=0; i<n; i++) for(int j=0; j<n; j++) ML_MAT(*U, i, j) = ML_MAT(A, i, j);

    for (int i = 0; i < n; i++) {
        // Find pivot
        int max_row = i;
        double max_val = fabs(ML_MAT(*U, i, i));
        for (int k = i + 1; k < n; k++) {
            if (fabs(ML_MAT(*U, k, i)) > max_val) {
                max_val = fabs(ML_MAT(*U, k, i));
                max_row = k;
            }
        }
        if (max_val == 0.0) return -1; // Singular

        // Swap rows in U, L, and P
        if (max_row != i) {
            for (int k = 0; k < n; k++) {
                double tmp = ML_MAT(*U, i, k); ML_MAT(*U, i, k) = ML_MAT(*U, max_row, k); ML_MAT(*U, max_row, k) = tmp;
                tmp = ML_MAT(*L, i, k); ML_MAT(*L, i, k) = ML_MAT(*L, max_row, k); ML_MAT(*L, max_row, k) = tmp;
            }
            int tmp = P[i]; P[i] = P[max_row]; P[max_row] = tmp;
        }

        // Elimination
        for (int k = i + 1; k < n; k++) {
            ML_MAT(*L, k, i) = ML_MAT(*U, k, i) / ML_MAT(*U, i, i);
            for (int j = i; j < n; j++) {
                ML_MAT(*U, k, j) -= ML_MAT(*L, k, i) * ML_MAT(*U, i, j);
            }
        }
    }
    return 0;
}

// Solve Ax = b using LU with Pivoting
inline int ml_matrix_solve_lu(ml_matrix A, double *b, double *x) {
    int n = A.rows;
    if (n <= 0) return -1;
    ml_matrix L, U;
    int *P = (int*)malloc(n * sizeof(int));
    if (ml_matrix_lu(A, &L, &U, P) != 0) { free(P); return -1; }

    // Permute b: Pb
    double *Pb = (double*)malloc(n * sizeof(double));
    for(int i=0; i<n; i++) Pb[i] = b[P[i]];

    // Forward substitution: Ly = Pb
    double *y = (double*)calloc((size_t)n, sizeof(double));
    for (int i = 0; i < n; i++) {
        double sum = 0;
        for (int j = 0; j < i; j++) sum += ML_MAT(L, i, j) * y[j];
        y[i] = Pb[i] - sum;
    }

    // Backward substitution: Ux = y
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0;
        for (int j = i + 1; j < n; j++) sum += ML_MAT(U, i, j) * x[j];
        x[i] = (y[i] - sum) / ML_MAT(U, i, i);
    }

    ml_matrix_free(L); ml_matrix_free(U); free(P); free(Pb); free(y);
    return 0;
}

// QR Decomposition (Householder Reflections)
// Returns 0 on success.
inline int ml_matrix_qr(ml_matrix A, ml_matrix *Q, ml_matrix *R) {
    int m = A.rows;
    int n = A.cols;
    *R = ml_matrix_create(m, n);
    *Q = ml_matrix_identity(m);

    // Copy A to R
    for(int i=0; i<m; i++) for(int j=0; j<n; j++) ML_MAT(*R, i, j) = ML_MAT(A, i, j);

    for (int k = 0; k < n; k++) {
        // Calculate Householder vector
        double norm_x = 0;
        for (int i = k; i < m; i++) norm_x += ML_MAT(*R, i, k) * ML_MAT(*R, i, k);
        norm_x = sqrt(norm_x);

        if (ML_MAT(*R, k, k) > 0) norm_x = -norm_x;

        double *v = (double*)calloc(m, sizeof(double));
        for (int i = k; i < m; i++) v[i] = ML_MAT(*R, i, k);
        v[k] -= norm_x;

        double norm_v = 0;
        for (int i = k; i < m; i++) norm_v += v[i] * v[i];
        if (norm_v == 0) { free(v); continue; }

        // Apply reflection to R from the LEFT: R = H_k * R
        for (int j = k; j < n; j++) {
            double dot = 0;
            for (int i = k; i < m; i++) dot += v[i] * ML_MAT(*R, i, j);
            for (int i = k; i < m; i++) ML_MAT(*R, i, j) -= 2 * v[i] * dot / norm_v;
        }

        // Apply reflection to Q from the RIGHT: Q = Q * H_k
        for (int i = 0; i < m; i++) {
            double dot = 0;
            for (int j = k; j < m; j++) dot += ML_MAT(*Q, i, j) * v[j];
            for (int j = k; j < m; j++) ML_MAT(*Q, i, j) -= 2 * dot * v[j] / norm_v;
        }
        free(v);
    }
    return 0;
}

// Power Iteration for Dominant Eigenvalue & Eigenvector
inline double ml_matrix_eigen_power(ml_matrix A, double *eigenvector, int max_iter) {
    int n = A.rows;
    double *b = (double*)calloc(n, sizeof(double));
    double *b_next = (double*)calloc(n, sizeof(double));
    for(int i=0; i<n; i++) b[i] = 1.0; // Initial guess

    double lambda = 0;
    for(int iter=0; iter<max_iter; iter++) {
        for(int i=0; i<n; i++) {
            b_next[i] = 0;
            for(int j=0; j<n; j++) b_next[i] += ML_MAT(A, i, j) * b[j];
        }
        double max_val = 0;
        for(int i=0; i<n; i++) if(fabs(b_next[i]) > fabs(max_val)) max_val = b_next[i];
        if(max_val == 0) break;
        for(int i=0; i<n; i++) b_next[i] /= max_val; // Normalize
        lambda = max_val;
        memcpy(b, b_next, n * sizeof(double));
    }
    memcpy(eigenvector, b, n * sizeof(double));
    free(b); free(b_next);
    return lambda;
}

// QR Algorithm with Wilkinson Shift for all Eigenvalues (Symmetric Matrices)
inline double* ml_matrix_eigen_qr(ml_matrix A, int max_iter) {
    int n = A.rows;
    ml_matrix Ak = ml_matrix_create(n, n);
    for(int i=0; i<n; i++) for(int j=0; j<n; j++) ML_MAT(Ak, i, j) = ML_MAT(A, i, j);

    for(int iter=0; iter<max_iter; iter++) {
        // True Wilkinson Shift (from bottom 2x2 principal submatrix)
        double a_w = ML_MAT(Ak, n-2, n-2);
        double b_w = ML_MAT(Ak, n-2, n-1);
        double c_w = ML_MAT(Ak, n-1, n-1);
        double delta_w = (a_w - c_w) / 2.0;
        double sign_delta = (delta_w >= 0) ? 1.0 : -1.0;
        double denom = delta_w + sign_delta * sqrt(delta_w * delta_w + b_w * b_w);
        double mu = (denom == 0.0) ? c_w : c_w - (b_w * b_w) / denom;
        ml_matrix I = ml_matrix_identity(n);
        for(int i=0; i<n; i++) ML_MAT(I, i, i) -= mu; // I = I - mu*I -> actually we need A - mu*I

        // Shift: Ak = Ak - mu*I
        for(int i=0; i<n; i++) ML_MAT(Ak, i, i) -= mu;

        ml_matrix Q, R;
        ml_matrix_qr(Ak, &Q, &R);
        ml_matrix Ak_next = ml_matrix_mul(R, Q);

        // Unshift: Ak = Ak_next + mu*I
        for(int i=0; i<n; i++) ML_MAT(Ak_next, i, i) += mu;

        ml_matrix_free(Ak);
        Ak = Ak_next;
        ml_matrix_free(Q); ml_matrix_free(R); ml_matrix_free(I);
    }

    double* eigenvalues = (double*)malloc(n * sizeof(double));
    for(int i=0; i<n; i++) eigenvalues[i] = ML_MAT(Ak, i, i);
    ml_matrix_free(Ak);
    return eigenvalues;
}

#endif
