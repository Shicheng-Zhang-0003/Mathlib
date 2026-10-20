#ifndef MATHLIB_ML_COMPLEXAN_H
#define MATHLIB_ML_COMPLEXAN_H

#include "ml_compiler.h"
#include "ml_core.h"
#include "ml_types.h"
#include "ml_complex.h"
#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * v12R2 COMPLEX ANALYSIS PACK
 *
 * Polynomial evaluation / all-roots (Durand-Kerner), complex elementary
 * functions beyond exp/log, numerical derivatives, contour and residue
 * calculus (Cauchy circle integral), complex quadratic roots and the
 * Mandelbrot escape-time iteration.
 * ========================================================================== */

typedef cplx (*ml_cplx_func_t) (cplx z, void *ctx);

/* c[0] is the coefficient of z^n (degree n), c[n] the constant term. */
ML_API cplx ml_cplx_poly_eval (const cplx *c, int n, cplx z);
ML_API ml_status_t ml_cplx_poly_roots (const cplx *c, int n, cplx *roots);

ML_API cplx ml_cplx_sin (cplx z);
ML_API cplx ml_cplx_cos (cplx z);
ML_API cplx ml_cplx_tan (cplx z);

ML_API cplx ml_cplx_derivative (ml_cplx_func_t f, cplx z, void *ctx);

/* Counter-clockwise contour integral of f along |z-center| = radius,
 * approximated by n trapezoidal nodes (spectrally accurate for analytic f). */
ML_API cplx ml_cplx_circle_integral (ml_cplx_func_t f, cplx center,
                                     double radius, int n, void *ctx);

/* Residue of f at an isolated pole, via the circle integral. */
ML_API cplx ml_cplx_residue (ml_cplx_func_t f, cplx pole, double radius, int n,
                             void *ctx);

/* Roots of a z^2 + b z + c = 0. */
ML_API void ml_cplx_solve_quadratic (cplx a, cplx b, cplx c, cplx *r1, cplx *r2);

/* Mandelbrot iteration count until |z| > radius (or max_iter reached). */
ML_API int ml_cplx_escape (cplx c, int max_iter, double radius);

#ifdef __cplusplus
}
#endif
#endif /* MATHLIB_ML_COMPLEXAN_H */
