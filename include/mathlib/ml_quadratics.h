#ifndef MATHLIB_ML_QUADRATICS_H
#define MATHLIB_ML_QUADRATICS_H

#include "ml_compiler.h"
#include "ml_core.h"
#ifdef __cplusplus
extern "C" {
#endif

ML_API double ml_equation(double a, double b, double c, double x);
ML_API double ml_formula_pos(double a, double b, double c);
ML_API double ml_formula_neg(double a, double b, double c);
/* DESPOT-AUDIT: returns number of distinct real roots written: 0 (none /
 * degenerate 0==0), 1 (one real, cubic casus single), 2 (double root:
 * quadratic Delta==0 or cubic Delta==0), 3 (three reals, casus
 * irreducibilis trig). r1/r2 are NaN when unused by the cubic path;
 * quadratic-delegation path leaves r2 untouched (legacy). */
ML_API int ml_cubic(double a, double b, double c, double d, double *r0, double *r1, double *r2);

#ifdef __cplusplus
}
#endif
#endif /* MATHLIB_ML_QUADRATICS_H */
