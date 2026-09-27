#ifndef MATHLIB_ML_ANALYTIC_NT_H
#define MATHLIB_ML_ANALYTIC_NT_H
#include "ml_compiler.h"
#include "ml_core.h"
#include "ml_complex.h"
/* ml_hurwitz_zeta: s>1 direct Euler-Maclaurin; a==1 delegates to zeta.
 * STUB: s<=1 with a!=1 returns NaN (analytic continuation not yet
 * implemented). */
ML_API cplx ml_hurwitz_zeta(double s, double a);
ML_API cplx ml_dirichlet_eta_cplx(cplx s);
ML_API double ml_theta3(double q);
ML_API double ml_partition_p(int n);
ML_API cplx ml_zeta_cplx(cplx s);
#endif
