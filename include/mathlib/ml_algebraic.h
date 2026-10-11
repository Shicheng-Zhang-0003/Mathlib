#ifndef MATHLIB_ML_ALGEBRAIC_H
#define MATHLIB_ML_ALGEBRAIC_H

#include <stdint.h>
#include "ml_compiler.h"
#include "ml_core.h"
#include "ml_types.h"
#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Computational algebra: Sylvester resultants and discriminants (exact over
 * int64 via the staged Bareiss determinant), Sturm real-root isolation
 * (long-double tier, degrees <= 8), Descartes' bound (exact), Specht
 * dimensions and S_n class sizes (hook formula, exact), Rodrigues so(3).
 *
 * Integer paths are exact-or-fail-loud (OVERFLOW / INVALID_ARG); Sturm is
 * LD with an exact-zero snap for small integer coefficients (documented
 * regression-guard tier, not a certified solver). No oracle ULP claimed.
 * ========================================================================== */

/* Elimination (coeffs ascending: c[0] + c[1] x + ...; |c| <= 1e6, deg <= 12) */
ML_API ml_status_t ml_resultant(const long long *p, int degp,
                                const long long *q, int degq, long long *res);
ML_API ml_status_t ml_discriminant(const long long *p, int deg, long long *disc);
/* Sturm: distinct real roots in (a, b]; -1 on invalid/degenerate. */
ML_API int ml_sturm_count(const long long *p, int deg, double a, double b);
/* Descartes: sign changes = positive-root upper bound (exact). */
ML_API int ml_descartes_bound(const long long *p, int deg);

/* Symmetric group (partition of n as Ferrers rows, n <= 20 for hook dims) */
ML_API uint64_t ml_hook_dimension(const int *parts, int len);
ML_API uint64_t ml_sn_class_size(const int *parts, int len);

/* so(3): R = I + sin t K + (1-cos t) K^2 for unit axis (validated). */
ML_API ml_status_t ml_rodrigues_exp(const double axis[3], double theta, double R[9]);

#ifdef __cplusplus
}
#endif
#endif /* MATHLIB_ML_ALGEBRAIC_H */
