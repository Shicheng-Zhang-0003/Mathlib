#ifndef MATHLIB_ML_MEASURE_H
#define MATHLIB_ML_MEASURE_H

#include "ml_compiler.h"
#include "ml_core.h"
#include "ml_types.h"
#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Measure-theoretic analysis: Lebesgue integral of simple functions and the
 * Dirichlet showcase ( rationals have measure zero), discrete L^p geometry,
 * and Hölder / Minkowski / Markov instance verifiers.
 *
 * Sample-based by design (no function closures in C99 without shared
 * state): f/g are uniform samples with spacing dx, so sums stand in for
 * integrals. Exactness: simple sums run in long double; verdicts carry a
 * caller tol (default 1e-9) and fail closed to 0. No oracle ULP claimed.
 * ========================================================================== */

ML_API double ml_lebesgue_dirichlet (double a, double b);
ML_API double ml_simple_integral (const double *values, const double *measures, int n);
ML_API double ml_lp_norm_data (const double *f, int n, double dx, double p);
ML_API double ml_seq_lp_norm (const double *a, int n, double p);
ML_API int ml_holder_check (const double *f, const double *g, int n,
                            double dx, double p, double tol);
ML_API int ml_minkowski_check (const double *f, const double *g, int n,
                               double dx, double p, double tol);
ML_API double ml_markov_bound (double mean_abs, double t);

#ifdef __cplusplus
}
#endif
#endif /* MATHLIB_ML_MEASURE_H */
