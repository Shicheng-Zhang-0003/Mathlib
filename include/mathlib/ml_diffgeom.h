#ifndef MATHLIB_ML_DIFFGEOM_H
#define MATHLIB_ML_DIFFGEOM_H

#include "ml_compiler.h"
#include "ml_core.h"
#include "ml_types.h"
#include "ml_numerical.h"
#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Differential geometry of curves and parametric surfaces, derivatives by
 * scale-aware central differences (same philosophy as ml_derivative).
 *
 * - Graph curvature reuses the staged ml_derivative/ml_second_derivative.
 * - Parametric Frenet frames need r', r'', r''' (5-point LD stencils).
 * - Surfaces take ml_surf_func_t callbacks; E,F,G / L,M,N feed Brioschi.
 * - Degenerate frames (straight motion, umbilic-free failure, r_u x r_v = 0)
 *   are SINGULAR with NaN-poisoned outs, never silent garbage.
 * No oracle ULP claimed.
 * ========================================================================== */

typedef double (*ml_surf_func_t) (double u, double v);

ML_API double ml_curvature_graph (ml_func_t f, double x);
ML_API double ml_curvature_param (ml_func_t fx, ml_func_t fy, double t);
ML_API ml_status_t ml_frenet_frame (ml_func_t fx, ml_func_t fy, ml_func_t fz,
                                    double t, double T [3], double N [3],
                                    double B [3], double *kappa, double *tau);
ML_API ml_status_t ml_surf_curvatures (ml_surf_func_t rx, ml_surf_func_t ry,
                                       ml_surf_func_t rz, double u, double v,
                                       double *K, double *H);

#ifdef __cplusplus
}
#endif
#endif /* MATHLIB_ML_DIFFGEOM_H */
