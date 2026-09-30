#ifndef MATHLIB_ML_KELVIN_H
#define MATHLIB_ML_KELVIN_H

#include "ml_core.h"
#include "ml_compiler.h"
#include "ml_integral.h"
#ifdef __cplusplus
extern "C" {
#endif

ML_API double ml_kelvin_ber(double x);
ML_API double ml_kelvin_bei(double x);
ML_API double ml_kelvin_ker(double x);
ML_API double ml_kelvin_kei(double x);

#ifdef __cplusplus
}
#endif
#endif
