#ifndef MATHLIB_ML_TOPOLOGY_H
#define MATHLIB_ML_TOPOLOGY_H

#include <stdint.h>
#include "ml_compiler.h"
#include "ml_core.h"
#include "ml_types.h"
#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Simplicial homology (mod 2) and exact integer linear algebra.
 *
 * Complexes are passed closed (every face present) as flattened vertex
 * lists: verts holds all vertex tuples back to back, off[k]..off[k+1] is
 * simplex k (strictly ascending vertices), sdim[k] its dimension.
 * Caps: vertices < 64, simplices <= 1024, dim <= 8. Non-closed or
 * malformed input -> INVALID_ARG (never a wrong Betti number).
 *
 * Bareiss fraction-free elimination is exact over int64 with entry cap
 * |a| <= 1e6 and dims <= 24; overflow past int64 reports OVERFLOW (-1 for
 * the rank query). No oracle ULP claimed.
 * ========================================================================== */

ML_API long long ml_euler_char (const int *counts, int dim);
ML_API ml_status_t ml_simplicial_betti_f2 (int nverts, const int *verts,
                                           const int *off, const int *sdim,
                                           int nsimps, int *betti, int maxdim);
ML_API int ml_chain_f2_is_complex (const int *mats, const int *rows,
                                   const int *cols, const int *offm, int nmaps);
ML_API ml_status_t ml_chain_betti_f2 (const int *mats, const int *rows,
                                      const int *cols, const int *offm,
                                      int nmaps, int *betti);
ML_API int ml_bareiss_rank (const long long *A, int m, int n);
ML_API ml_status_t ml_bareiss_det (const long long *A, int n, long long *det);

#ifdef __cplusplus
}
#endif
#endif /* MATHLIB_ML_TOPOLOGY_H */
