#ifndef MATHLIB_ML_QUATERNION_H
#define MATHLIB_ML_QUATERNION_H

#include "ml_compiler.h"
#include "ml_trig.h"
#include "ml_core.h"

typedef struct { double w, x, y, z; } ml_quat;

ML_API ml_quat ml_quat_mul(ml_quat a, ml_quat b);
ML_API ml_quat ml_quat_normalize(ml_quat q);
ML_API ml_quat ml_quat_slerp(ml_quat a, ml_quat b, double t);
ML_API ml_quat ml_quat_conjugate(ml_quat q);
ML_API ml_quat ml_quat_inverse(ml_quat q);
ML_API ml_quat ml_quat_from_axis_angle(double ax, double ay, double az, double angle);

#endif /* MATHLIB_ML_QUATERNION_H */
