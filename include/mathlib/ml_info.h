#ifndef MATHLIB_ML_INFO_H
#define MATHLIB_ML_INFO_H
#include "ml_compiler.h"
#include "ml_core.h"
/* Natural-log (nats) information measures. p/q are un-normalized
 * non-negative weights (no simplex renormalization is applied); entropy is
 * -sum p log p, KL is sum p log(p/q) with the 0*log0=0 convention. */
ML_API double ml_entropy(const double *p, int n);
ML_API double ml_kl_div(const double *p, const double *q, int n);
ML_API double ml_cross_entropy(const double *p, const double *q, int n);
/* Discrete mutual information over an nr x nc joint table (row-major).
 * 1 <= nr, nc <= 1024. Returns +Inf for zero-marginal support mismatch. */
ML_API double ml_mi_discrete(const double *joint, int nr, int nc);
ML_API double ml_logistic(double x);
ML_API double ml_softplus(double x);
#endif
