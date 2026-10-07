#ifndef MATHLIB_ML_ANALYTIC_NT_H
#define MATHLIB_ML_ANALYTIC_NT_H
#include "ml_compiler.h"
#include "ml_core.h"
#include "ml_complex.h"
#ifdef __cplusplus
extern "C" {
#endif
/* ml_hurwitz_zeta: s>0 direct LD Euler-Maclaurin (N=200, B2..B12 tail);
 * a==1 delegates to zeta. s<=0 with a!=1 remains a fail-loud NaN STUB
 * (needs the Hurwitz functional equation; documented in
 * KNOWN_LIMITATIONS.md). Previously the whole s<=1,a!=1 region was a
 * stub; Round-4 closed 0<s<=1 via the same EM analytic continuation. */
ML_API cplx ml_hurwitz_zeta(double s, double a);
ML_API cplx ml_dirichlet_eta_cplx(cplx s);
ML_API double ml_theta3(double q);
ML_API double ml_partition_p(int n);
ML_API cplx ml_zeta_cplx(cplx s);
#ifdef __cplusplus
}
#endif
#endif
