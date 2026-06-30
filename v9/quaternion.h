#ifndef MATHLIB_QUATERNION_H
#define MATHLIB_QUATERNION_H

#include "profiles.h"
#include "ml_core.h"

typedef struct { double w, x, y, z; } ml_quat;

static inline ml_quat ml_quat_mul(ml_quat a, ml_quat b) {
    return (ml_quat){
        a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z,
        a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y,
        a.w*b.y - a.x*b.z + a.y*b.w + a.z*b.x,
        a.w*b.z + a.x*b.y - a.y*b.x + a.z*b.w
    };
}

static inline ml_quat ml_quat_normalize(ml_quat q) {
    double mag2 = q.w*q.w + q.x*q.x + q.y*q.y + q.z*q.z;
    double inv_mag = ml_rsqrt(mag2);
    return (ml_quat){q.w*inv_mag, q.x*inv_mag, q.y*inv_mag, q.z*inv_mag};
}

static inline ml_quat ml_quat_slerp(ml_quat a, ml_quat b, double t) {
    double dot = a.w*b.w + a.x*b.x + a.y*b.y + a.z*b.z;
    if (dot < 0) { b.w=-b.w; b.x=-b.x; b.y=-b.y; b.z=-b.z; dot=-dot; }
    if (dot > 0.9995) {
        return ml_quat_normalize((ml_quat){
            a.w + t*(b.w-a.w), a.x + t*(b.x-a.x),
            a.y + t*(b.y-a.y), a.z + t*(b.z-a.z)
        });
    }
    double theta = acos(dot);
    double sin_theta = ml_sin(theta); // Profile-aware trig
    double s0 = ml_sin((1.0 - t) * theta) / sin_theta;
    double s1 = ml_sin(t * theta) / sin_theta;
    return (ml_quat){
        s0*a.w + s1*b.w, s0*a.x + s1*b.x,
        s0*a.y + s1*b.y, s0*a.z + s1*b.z
    };
}
#endif
