#ifndef MATHLIB_ML_FIXED_POINT_H
#define MATHLIB_ML_FIXED_POINT_H

#include <stdint.h>
#include "ml_compiler.h"

typedef int32_t ml_q16_16_t;

#define ML_FIXED_PI 205887
#define ML_FIXED_HALF_PI 102943
#define ML_FIXED_TWO_PI 411774
#define ML_FIXED_CORDIC_GAIN 39797

ML_API ml_q16_16_t ml_fixed_mul(ml_q16_16_t a, ml_q16_16_t b);
ML_API void ml_cordic_sincos_fixed(ml_q16_16_t theta, ml_q16_16_t *sin_out, ml_q16_16_t *cos_out);

#endif /* MATHLIB_ML_FIXED_POINT_H */
