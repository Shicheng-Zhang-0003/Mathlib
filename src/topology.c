#include "ml_compiler.h"
#include "ml_topology.h"
#include "ml_core.h"
#include <stddef.h>
#include <stdlib.h>

/* Simplicial homology over F2 + exact Bareiss integer rank/det. */

#define MLT_VMAX 64
#define MLT_SMAX 1024
#define MLT_DMAX 8
#define MLT_ENTRY_MAX 1000000LL

ML_API long long ml_euler_char(const int *counts, int dim) {
    long long chi = 0;
    int k;
    if (!counts || dim < 0 || dim > MLT_DMAX) return INT64_MIN;
    for (k = 0; k <= dim; k++) {
        long long c = (long long)counts[k];
        if (c < 0 || c > MLT_SMAX) return INT64_MIN;
        if ((k & 1) == 0) {
            if (chi > INT64_MAX - c) return INT64_MIN;
            chi += c;
        } else {
            if (chi < INT64_MIN + c) return INT64_MIN;
            chi -= c;
        }
    }
    return chi;
}

/* GF(2) row rank of a rows x cols 0/1 matrix (row-major ints).
 * Heap bitsets (a 1024x1024 bitset is 128KB: never stack). */
static int mlt_rank_f2(const int *mat, int rows, int cols) {
    uint64_t (*work)[16] = NULL;
    int r = 0;
    int c;
    int i;
    int j;
    if (rows <= 0 || cols <= 0 || cols > 1024 || rows > 1024) return -1;
    if ((size_t)rows > (size_t)-1 / (16 * sizeof(uint64_t))) return -1;
    work = (uint64_t (*)[16])malloc((size_t)rows * 16 * sizeof(uint64_t));
    if (!work) return -1;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < 16; j++) work[i][j] = 0;
        for (j = 0; j < cols; j++) {
            int v = mat[(size_t)i * (size_t)cols + (size_t)j];
            if (v != 0 && v != 1) return -1;
            if (v) work[i][(size_t)j / 64U] |= 1ULL << ((unsigned)j % 64U);
        }
    }
    c = 0;
    for (i = 0; i < rows && c < cols; i++) {
        int piv = -1;
        int k;
        for (k = i; k < rows; k++) {
            if ((work[k][(size_t)c / 64U] >> ((unsigned)c % 64U)) & 1ULL) {
                piv = k;
                break;
            }
        }
        if (piv < 0) {
            /* Find next set column in rows i.. (full pivot search). */
            int cc;
            int found = 0;
            for (cc = c + 1; cc < cols && !found; cc++) {
                for (k = i; k < rows; k++) {
                    if ((work[k][(size_t)cc / 64U] >> ((unsigned)cc % 64U)) & 1ULL) {
                        piv = k;
                        c = cc;
                        found = 1;
                        break;
                    }
                }
            }
            if (!found) break;
        }
        if (piv != i) {
            for (j = 0; j < 16; j++) {
                uint64_t t = work[i][j];
                work[i][j] = work[piv][j];
                work[piv][j] = t;
            }
        }
        for (k = i + 1; k < rows; k++) {
            if ((work[k][(size_t)c / 64U] >> ((unsigned)c % 64U)) & 1ULL) {
                for (j = 0; j < 16; j++) work[k][j] ^= work[i][j];
            }
        }
        r++;
        c++;
    }
    free(work);
    return r;
}

ML_API ml_status_t ml_simplicial_betti_f2(int nverts, const int *verts,
                                          const int *off, const int *sdim,
                                          int nsimps, int *betti, int maxdim) {
    /* Validate, bucket by dim, check closure, build boundary maps, rank. */
    int counts[9];
    int *idx_of = NULL;
    int k;
    int s;
    if (nverts < 1 || nverts > MLT_VMAX) return ML_ERR_INVALID_ARG;
    if (!verts || !off || !sdim || !betti || nsimps < 1 || nsimps > MLT_SMAX)
        return ML_ERR_INVALID_ARG;
    if (maxdim < 0 || maxdim > MLT_DMAX) return ML_ERR_INVALID_ARG;
    for (k = 0; k < 9; k++) counts[k] = 0;
    /* Validate simplices: ascending verts in range, dim matches length. */
    for (s = 0; s < nsimps; s++) {
        int len;
        int j;
        if (sdim[s] < 0 || sdim[s] > maxdim) return ML_ERR_INVALID_ARG;
        if (off[s] < 0 || off[s + 1] < off[s]) return ML_ERR_INVALID_ARG;
        len = off[s + 1] - off[s];
        if (len != sdim[s] + 1) return ML_ERR_INVALID_ARG;
        for (j = 0; j < len; j++) {
            int v = verts[off[s] + j];
            if (v < 0 || v >= nverts) return ML_ERR_INVALID_ARG;
            if (j > 0 && v <= verts[off[s] + j - 1]) return ML_ERR_INVALID_ARG;
        }
        counts[sdim[s]]++;
    }
    /* Index simplices by (dim, position) for boundary construction. */
    if ((size_t)nsimps > (size_t)-1 / sizeof(int)) return ML_ERR_INVALID_ARG;
    idx_of = (int *)malloc((size_t)nsimps * sizeof(int));
    if (!idx_of) return ML_ERR_WORKSPACE;
    {
        int pos[9];
        for (k = 0; k < 9; k++) pos[k] = 0;
        for (s = 0; s < nsimps; s++) idx_of[s] = pos[sdim[s]]++;
    }
    /* Closure: every face of every simplex must exist. */
    for (s = 0; s < nsimps; s++) {
        int d = sdim[s];
        int f;
        if (d == 0) continue;
        for (f = 0; f <= d; f++) {
            int t;
            int found = 0;
            for (t = 0; t < nsimps; t++) {
                int td = sdim[t];
                int j2;
                int match;
                if (td != d - 1) continue;
                if (off[t + 1] - off[t] != d) continue;
                match = 1;
                for (j2 = 0; j2 < d; j2++) {
                    int want = verts[off[s] + j2 + (j2 >= f ? 1 : 0)];
                    if (verts[off[t] + j2] != want) { match = 0; break; }
                }
                if (match) { found = 1; break; }
            }
            if (!found) {
                free(idx_of);
                return ML_ERR_INVALID_ARG;
            }
        }
    }
    /* Betti_k = C_k - rank(d_k) - rank(d_{k+1}). */
    for (k = 0; k <= maxdim; k++) {
        int rk = 0;
        int rk1 = 0;
        if (k > 0 && counts[k] > 0 && counts[k - 1] > 0) {
            /* Boundary d_k: rows = (k-1)-simplices, cols = k-simplices. */
            int rows = counts[k - 1];
            int cols = counts[k];
            int *mat = NULL;
            int t;
            if ((size_t)rows * (size_t)cols > (size_t)1048576) {
                free(idx_of);
                return ML_ERR_WORKSPACE;
            }
            mat = (int *)malloc((size_t)rows * (size_t)cols * sizeof(int));
            if (!mat) {
                free(idx_of);
                return ML_ERR_WORKSPACE;
            }
            for (t = 0; t < rows * cols; t++) mat[t] = 0;
            for (t = 0; t < nsimps; t++) {
                int d;
                int f;
                if (sdim[t] != k) continue;
                d = k;
                for (f = 0; f <= d; f++) {
                    int u;
                    for (u = 0; u < nsimps; u++) {
                        int j2;
                        int match;
                        if (sdim[u] != k - 1) continue;
                        match = 1;
                        for (j2 = 0; j2 < d; j2++) {
                            int want = verts[off[t] + j2 + (j2 >= f ? 1 : 0)];
                            if (verts[off[u] + j2] != want) { match = 0; break; }
                        }
                        if (match) {
                            mat[(size_t)idx_of[u] * (size_t)cols + (size_t)idx_of[t]] ^= 1;
                            break;
                        }
                    }
                }
            }
            rk = mlt_rank_f2(mat, rows, cols);
            free(mat);
            if (rk < 0) {
                free(idx_of);
                return ML_ERR_INTERNAL;
            }
        }
        if (k < maxdim && counts[k + 1] > 0 && counts[k] > 0) {
            int rows = counts[k];
            int cols = counts[k + 1];
            int *mat = NULL;
            int t;
            if ((size_t)rows * (size_t)cols > (size_t)1048576) {
                free(idx_of);
                return ML_ERR_WORKSPACE;
            }
            mat = (int *)malloc((size_t)rows * (size_t)cols * sizeof(int));
            if (!mat) {
                free(idx_of);
                return ML_ERR_WORKSPACE;
            }
            for (t = 0; t < rows * cols; t++) mat[t] = 0;
            for (t = 0; t < nsimps; t++) {
                int d;
                int f;
                if (sdim[t] != k + 1) continue;
                d = k + 1;
                for (f = 0; f <= d; f++) {
                    int u;
                    for (u = 0; u < nsimps; u++) {
                        int j2;
                        int match;
                        if (sdim[u] != k) continue;
                        match = 1;
                        for (j2 = 0; j2 <= d - 1; j2++) {
                            int want = verts[off[t] + j2 + (j2 >= f ? 1 : 0)];
                            if (verts[off[u] + j2] != want) { match = 0; break; }
                        }
                        if (match) {
                            mat[(size_t)idx_of[u] * (size_t)cols + (size_t)idx_of[t]] ^= 1;
                            break;
                        }
                    }
                }
            }
            rk1 = mlt_rank_f2(mat, rows, cols);
            free(mat);
            if (rk1 < 0) {
                free(idx_of);
                return ML_ERR_INTERNAL;
            }
        }
        betti[k] = counts[k] - rk - rk1;
        if (betti[k] < 0) {
            free(idx_of);
            return ML_ERR_INTERNAL;
        }
    }
    free(idx_of);
    return ML_SUCCESS;
}

ML_API int ml_chain_f2_is_complex(const int *mats, const int *rows,
                                  const int *cols, const int *offm, int nmaps) {
    /* d_{i} d_{i+1} = 0 over F2 for each adjacent pair. */
    int i;
    if (!mats || !rows || !cols || !offm || nmaps < 1 || nmaps > 8) return 0;
    for (i = 0; i + 1 < nmaps; i++) {
        int r1 = rows[i], c1 = cols[i];
        int r2 = rows[i + 1], c2 = cols[i + 1];
        int a;
        int b;
        if (r1 <= 0 || c1 <= 0 || r2 <= 0 || c2 <= 0) return 0;
        if (r1 > 1024 || c1 > 1024 || r2 > 1024 || c2 > 1024) return 0;
        if (c1 != r2) return 0;
        for (a = 0; a < r1; a++) {
            for (b = 0; b < c2; b++) {
                int s = 0;
                int k;
                for (k = 0; k < c1; k++) {
                    int x = mats[offm[i] + a * c1 + k];
                    int y = mats[offm[i + 1] + k * c2 + b];
                    if ((x != 0 && x != 1) || (y != 0 && y != 1)) return 0;
                    s ^= (x & y);
                }
                if (s != 0) return 0;
            }
        }
    }
    return 1;
}

ML_API ml_status_t ml_chain_betti_f2(const int *mats, const int *rows,
                                     const int *cols, const int *offm,
                                     int nmaps, int *betti) {
    /* H_k = ker(d_k)/im(d_{k+1}): dim C_k - rank(d_k) - rank(d_{k+1}).
     * Chain: C_n --d_n--> ... --d_1--> C_0, mats[i] is d_{n-i}... caller
     * passes maps in order d_1..d_n with dims: C_k has dim rows/cols shared
     * with neighbours (validated). betti has nmaps+1 entries. */
    int k;
    int *ranks = NULL;
    if (!mats || !rows || !cols || !offm || !betti || nmaps < 1 || nmaps > 8)
        return ML_ERR_INVALID_ARG;
    for (k = 0; k < nmaps; k++) {
        if (rows[k] <= 0 || cols[k] <= 0 || rows[k] > 1024 || cols[k] > 1024)
            return ML_ERR_INVALID_ARG;
        if (k > 0 && cols[k - 1] != rows[k]) return ML_ERR_INVALID_ARG;
    }
    if (!ml_chain_f2_is_complex(mats, rows, cols, offm, nmaps))
        return ML_ERR_INVALID_ARG;
    if ((size_t)nmaps > (size_t)-1 / sizeof(int)) return ML_ERR_INVALID_ARG;
    ranks = (int *)malloc((size_t)nmaps * sizeof(int));
    if (!ranks) return ML_ERR_WORKSPACE;
    for (k = 0; k < nmaps; k++) {
        ranks[k] = mlt_rank_f2(mats + offm[k], rows[k], cols[k]);
        if (ranks[k] < 0) {
            free(ranks);
            return ML_ERR_INTERNAL;
        }
    }
    /* C_0 dim = rows[0]; C_k dim = cols[k-1] = rows[k]. */
    for (k = 0; k <= nmaps; k++) {
        int dimc;
        int rout;
        int rin;
        if (k == 0) {
            dimc = rows[0];
            rout = ranks[0];
            rin = 0;
        } else if (k == nmaps) {
            dimc = cols[nmaps - 1];
            rout = 0;
            rin = ranks[nmaps - 1];
        } else {
            dimc = rows[k];
            rout = ranks[k];
            rin = ranks[k - 1];
        }
        betti[k] = dimc - rout - rin;
        if (betti[k] < 0) {
            free(ranks);
            return ML_ERR_INTERNAL;
        }
    }
    free(ranks);
    return ML_SUCCESS;
}

ML_API int ml_bareiss_rank(const long long *A, int m, int n) {
    /* Fraction-free Bareiss, exact over int64. -1 on invalid/overflow. */
    long long *b = NULL;
    long long prev = 1;
    int rank = 0;
    int k;
    int i;
    int j;
    if (!A || m <= 0 || n <= 0 || m > 24 || n > 24) return -1;
    for (i = 0; i < m * n; i++) {
        if (A[i] > MLT_ENTRY_MAX || A[i] < -MLT_ENTRY_MAX) return -1;
    }
    if ((size_t)m * (size_t)n > (size_t)-1 / sizeof(long long)) return -1;
    b = (long long *)malloc((size_t)m * (size_t)n * sizeof(long long));
    if (!b) return -1;
    for (i = 0; i < m * n; i++) b[i] = A[i];
    for (k = 0; k < m && k < n; k++) {
        int piv = -1;
        for (i = k; i < m; i++) {
            if (b[(size_t)i * (size_t)n + (size_t)k] != 0) { piv = i; break; }
        }
        if (piv < 0) continue;
        if (piv != k) {
            for (j = k; j < n; j++) {
                long long t = b[(size_t)k * (size_t)n + (size_t)j];
                b[(size_t)k * (size_t)n + (size_t)j] = b[(size_t)piv * (size_t)n + (size_t)j];
                b[(size_t)piv * (size_t)n + (size_t)j] = t;
            }
        }
        for (i = k + 1; i < m; i++) {
            for (j = k + 1; j < n; j++) {
                long long aik = b[(size_t)i * (size_t)n + (size_t)k];
                long double num;
                long double q;
                if (aik == 0) continue;
                /* Exact in long double below 2^64 (guarded); Bareiss
                 * division is provably even, checked here. */
                num = (long double)b[(size_t)i * (size_t)n + (size_t)j] *
                          (long double)b[(size_t)k * (size_t)n + (size_t)k] -
                      (long double)aik *
                          (long double)b[(size_t)k * (size_t)n + (size_t)j];
                if (num >= 18446744073709551616.0L || num <= -18446744073709551616.0L) {
                    free(b);
                    return -1;
                }
                if (prev == 0) {
                    free(b);
                    return -1;
                }
                q = num / (long double)prev;
                if (q > 9223372036854775807.0L || q < -9223372036854775808.0L) {
                    free(b);
                    return -1;
                }
                if (q != (long double)(long long)q) {
                    free(b);
                    return -1;
                }
                b[(size_t)i * (size_t)n + (size_t)j] = (long long)q;
            }
        }
        for (i = k + 1; i < m; i++) b[(size_t)i * (size_t)n + (size_t)k] = 0;
        prev = b[(size_t)k * (size_t)n + (size_t)k];
        if (prev == 0) {
            free(b);
            return -1;
        }
        rank++;
    }
    free(b);
    return rank;
}

ML_API ml_status_t ml_bareiss_det(const long long *A, int n, long long *det) {
    long long *b = NULL;
    long long prev = 1;
    long long sign = 1;
    int k;
    int i;
    int j;
    if (ML_UNLIKELY(!A || !det)) return ML_ERR_INVALID_ARG;
    if (n <= 0 || n > 24) return ML_ERR_INVALID_ARG;
    for (i = 0; i < n * n; i++) {
        if (A[i] > MLT_ENTRY_MAX || A[i] < -MLT_ENTRY_MAX) return ML_ERR_INVALID_ARG;
    }
    if ((size_t)n * (size_t)n > (size_t)-1 / sizeof(long long)) return ML_ERR_INVALID_ARG;
    b = (long long *)malloc((size_t)n * (size_t)n * sizeof(long long));
    if (!b) return ML_ERR_WORKSPACE;
    for (i = 0; i < n * n; i++) b[i] = A[i];
    for (k = 0; k < n - 1; k++) {
        int piv = -1;
        for (i = k; i < n; i++) {
            if (b[(size_t)i * (size_t)n + (size_t)k] != 0) { piv = i; break; }
        }
        if (piv < 0) {
            *det = 0;
            free(b);
            return ML_SUCCESS;
        }
        if (piv != k) {
            for (j = k; j < n; j++) {
                long long t = b[(size_t)k * (size_t)n + (size_t)j];
                b[(size_t)k * (size_t)n + (size_t)j] = b[(size_t)piv * (size_t)n + (size_t)j];
                b[(size_t)piv * (size_t)n + (size_t)j] = t;
            }
            sign = -sign;
        }
        for (i = k + 1; i < n; i++) {
            for (j = k + 1; j < n; j++) {
                long long aik = b[(size_t)i * (size_t)n + (size_t)k];
                long long akj = b[(size_t)k * (size_t)n + (size_t)j];
                long long akk = b[(size_t)k * (size_t)n + (size_t)k];
                long long aij = b[(size_t)i * (size_t)n + (size_t)j];
                long double num = (long double)aij * (long double)akk -
                                  (long double)aik * (long double)akj;
                long double q;
                if (num >= 18446744073709551616.0L || num <= -18446744073709551616.0L) {
                    free(b);
                    return ML_ERR_OVERFLOW;
                }
                if (prev == 0) {
                    free(b);
                    return ML_ERR_INTERNAL;
                }
                q = num / (long double)prev;
                if (q > 9223372036854775807.0L || q < -9223372036854775808.0L) {
                    free(b);
                    return ML_ERR_OVERFLOW;
                }
                if (q != (long double)(long long)q) {
                    free(b);
                    return ML_ERR_INTERNAL;
                }
                b[(size_t)i * (size_t)n + (size_t)j] = (long long)q;
            }
        }
        prev = b[(size_t)k * (size_t)n + (size_t)k];
        if (prev == 0) {
            free(b);
            return ML_ERR_INTERNAL;
        }
    }
    *det = sign * b[(size_t)(n - 1) * (size_t)n + (size_t)(n - 1)];
    free(b);
    return ML_SUCCESS;
}
