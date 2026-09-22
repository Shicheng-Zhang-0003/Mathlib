#ifndef MATHLIB_ML_CALCULUS_H
#define MATHLIB_ML_CALCULUS_H

#include "ml_compiler.h"
#include "ml_core.h"
#include "ml_types.h"
#include "ml_numerical.h"

/* ============================================================================
 * Undergraduate calculus: quadrature theorems, interpolation, vector calc.
 *
 * Theorems: Gauss-Legendre exactness deg 2n-1, Lagrange existence/uniqueness,
 * cubic spline C2 minimal-curvature, FTC-backed adaptive quadrature,
 * Green/Stokes/Divergence via discrete div/curl consistency.
 *
 * Accuracy: nodes/weights in long double then rounded (<1 ULP vs mpmath
 * for n<=32). Interpolants use long-double barycentric accumulation.
 * ========================================================================== */

ML_API ml_status_t ml_gauss_legendre(int n, double *xs, double *ws);
ML_API double ml_gauss_legendre_integral(ml_func_t f, double a, double b, int n);
ML_API double ml_lagrange_interp(const double *x, const double *y, int n, double x0);
ML_API ml_status_t ml_cubic_spline_natural(const double *x, const double *y, int n,
                                           double *m2);
ML_API double ml_cubic_spline_eval(const double *x, const double *y,
                                   const double *m2, int n, double x0);
ML_API void ml_grad3(ml_func_t f3[3], double x, double y, double z, double h,
                     double *gx, double *gy, double *gz);
ML_API double ml_div3(ml_func_t f3[3], double x, double y, double z, double h);
ML_API void ml_curl3(ml_func_t f3[3], double x, double y, double z, double h,
                     double *cx, double *cy, double *cz);

#endif /* MATHLIB_ML_CALCULUS_H */
