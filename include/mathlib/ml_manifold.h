#ifndef MATHLIB_ML_MANIFOLD_H
#define MATHLIB_ML_MANIFOLD_H
#include "ml_compiler.h"
#include "ml_core.h"
#include "ml_types.h"
ML_API ml_status_t ml_proj_sphere(const double *x, double *y, int n);
/* Thin-QR Stiefel projection (MGS + one re-orthogonalization pass).
 * Orthonormalizes A's columns; NOT the polar-factor nearest point in
 * Frobenius norm. Requires full column rank: rank-deficient A returns
 * ML_ERR_SINGULAR. m >= n, m <= 32, n <= 32. */
ML_API ml_status_t ml_proj_stiefel(const double *A, double *Q, int m, int n);
ML_API ml_status_t ml_exp_sphere(const double *x, const double *v, double *y, int n);
ML_API double ml_sphere_dist(const double *x, const double *y, int n);
#endif
