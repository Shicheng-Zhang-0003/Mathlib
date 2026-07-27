#ifndef LIBMATHC_ML_COMPLEX_H
#define LIBMATHC_ML_COMPLEX_H

#include "ml_compiler.h"
#include "ml_core.h"
#include "ml_exp_log.h"
#include "ml_trig.h"

typedef struct { double real; double imag; } cplx;

/* Trivial inlines kept in header */
ML_INLINE cplx ml_cplx_add(cplx a, cplx b) { return (cplx){a.real + b.real, a.imag + b.imag}; }
ML_INLINE cplx ml_cplx_sub(cplx a, cplx b) { return (cplx){a.real - b.real, a.imag - b.imag}; }
ML_INLINE cplx ml_cplx_mul(cplx a, cplx b) { return (cplx){a.real * b.real - a.imag * b.imag, a.real * b.imag + a.imag * b.real}; }
ML_INLINE cplx ml_cplx_conjugate(cplx a) { return (cplx){a.real, -a.imag}; }

/* Exported complex operations */
ML_API double ml_cplx_abs(cplx a);
ML_API double ml_cplx_arg(cplx a);
ML_API cplx ml_cplx_div(cplx a, cplx b);
ML_API cplx ml_cplx_exponential(cplx a);
ML_API cplx ml_cplx_logarithm(cplx a);
ML_API cplx ml_cplx_power(cplx a, cplx b);

#endif /* LIBMATHC_ML_COMPLEX_H */
