#ifndef MATHLIB_ML_FIXED_POINT_H
#define MATHLIB_ML_FIXED_POINT_H

#include <stdint.h>
#include "ml_compiler.h"
#ifdef __cplusplus
extern "C" {
#endif

typedef int32_t ml_q16_16_t;

#define ML_FIXED_PI 205887
/* DESPOT-AUDIT: HALF_PI/TWO_PI were truncated, not rounded (off by 1 LSB
 * = 1.5e-5 rad bias). Correct: round(pi/2*65536)=102944,
 * round(2pi*65536)=411775. PI*65536=205887.41 truncates correctly. */
#define ML_FIXED_HALF_PI 102944
#define ML_FIXED_TWO_PI 411775
#define ML_FIXED_CORDIC_GAIN 39797

ML_API ml_q16_16_t ml_fixed_mul(ml_q16_16_t a, ml_q16_16_t b);
ML_API void ml_cordic_sincos_fixed(ml_q16_16_t theta, ml_q16_16_t *sin_out, ml_q16_16_t *cos_out);

#ifdef __cplusplus
}
#endif
#endif /* MATHLIB_ML_FIXED_POINT_H */
