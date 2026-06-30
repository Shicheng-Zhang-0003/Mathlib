#ifndef LIBMATHC_LINEAR_ALGEBRA_H
#define LIBMATHC_LINEAR_ALGEBRA_H

#include "ml_core.h"
#include <stdio.h>

// --- Struct Definitions ---
typedef struct { double x, y; } vec2;
typedef struct { double x, y, z; } vec3;
typedef struct { double m[9]; } mat3x3;

// --- Vec2 Operations ---
static inline vec2 vec2_add(vec2 a, vec2 b) { return (vec2){a.x + b.x, a.y + b.y}; }
static inline vec2 vec2_sub(vec2 a, vec2 b) { return (vec2){a.x - b.x, a.y - b.y}; }
static inline vec2 vec2_scale(vec2 a, double s) { return (vec2){a.x * s, a.y * s}; }
static inline double vec2_dot(vec2 a, vec2 b) { return a.x * b.x + a.y * b.y; }
static inline double vec2_cross(vec2 a, vec2 b) { return a.x * b.y - a.y * b.x; }
static inline double vec2_mag(vec2 a) { return ml_sqrt(a.x * a.x + a.y * a.y); }
static inline vec2 vec2_normalize(vec2 a) {
    double m = vec2_mag(a);
    if (m == 0.0) return (vec2){0.0, 0.0};
    return (vec2){a.x / m, a.y / m};
}

// --- Vec3 Operations ---
static inline vec3 vec3_add(vec3 a, vec3 b) { return (vec3){a.x + b.x, a.y + b.y, a.z + b.z}; }
static inline vec3 vec3_sub(vec3 a, vec3 b) { return (vec3){a.x - b.x, a.y - b.y, a.z - b.z}; }
static inline vec3 vec3_scale(vec3 a, double s) { return (vec3){a.x * s, a.y * s, a.z * s}; }
static inline double vec3_dot(vec3 a, vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static inline vec3 vec3_cross(vec3 a, vec3 b) {
    return (vec3){
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
}
static inline double vec3_mag(vec3 a) { return ml_sqrt(a.x * a.x + a.y * a.y + a.z * a.z); }
static inline vec3 vec3_normalize(vec3 a) {
    double m = vec3_mag(a);
    if (m == 0.0) return (vec3){0.0, 0.0, 0.0};
    return (vec3){a.x / m, a.y / m, a.z / m};
}

// --- Mat3x3 Operations ---
static inline mat3x3 mat3x3_identity() {
    return (mat3x3){{1,0,0, 0,1,0, 0,0,1}};
}

static inline mat3x3 mat3x3_mul(mat3x3 a, mat3x3 b) {
    mat3x3 out;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            out.m[i * 3 + j] = 0;
            for (int k = 0; k < 3; k++) {
                out.m[i * 3 + j] += a.m[i * 3 + k] * b.m[k * 3 + j];
            }
        }
    }
    return out;
}

static inline vec3 mat3x3_mul_vec3(mat3x3 m, vec3 v) {
    return (vec3){
        m.m[0]*v.x + m.m[1]*v.y + m.m[2]*v.z,
        m.m[3]*v.x + m.m[4]*v.y + m.m[5]*v.z,
        m.m[6]*v.x + m.m[7]*v.y + m.m[8]*v.z
    };
}

static inline double mat3x3_det(mat3x3 m) {
    return m.m[0] * (m.m[4] * m.m[8] - m.m[5] * m.m[7])
         - m.m[1] * (m.m[3] * m.m[8] - m.m[5] * m.m[6])
         + m.m[2] * (m.m[3] * m.m[7] - m.m[4] * m.m[6]);
}

static inline mat3x3 mat3x3_transpose(mat3x3 m) {
    mat3x3 out;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            out.m[j * 3 + i] = m.m[i * 3 + j];
        }
    }
    return out;
}

static inline mat3x3 mat3x3_inverse(mat3x3 m) {
    double det = mat3x3_det(m);
    if (ml_fabs(det) < 1e-12) {
        mat3x3 nan_mat;
        for(int i=0; i<9; i++) nan_mat.m[i] = 0.0/0.0;
        return nan_mat;
    }
    double inv_det = 1.0 / det;
    mat3x3 out;
    out.m[0] =  (m.m[4] * m.m[8] - m.m[5] * m.m[7]) * inv_det;
    out.m[1] = -(m.m[1] * m.m[8] - m.m[2] * m.m[7]) * inv_det;
    out.m[2] =  (m.m[1] * m.m[5] - m.m[2] * m.m[4]) * inv_det;
    out.m[3] = -(m.m[3] * m.m[8] - m.m[5] * m.m[6]) * inv_det;
    out.m[4] =  (m.m[0] * m.m[8] - m.m[2] * m.m[6]) * inv_det;
    out.m[5] = -(m.m[0] * m.m[5] - m.m[2] * m.m[3]) * inv_det;
    out.m[6] =  (m.m[3] * m.m[7] - m.m[4] * m.m[6]) * inv_det;
    out.m[7] = -(m.m[0] * m.m[7] - m.m[1] * m.m[6]) * inv_det;
    out.m[8] =  (m.m[0] * m.m[4] - m.m[1] * m.m[3]) * inv_det;
    return out;
}

static inline vec3 linear_solve_3x3(mat3x3 m, vec3 v) {
    mat3x3 inv = mat3x3_inverse(m);
    if (ml_isnan(inv.m[0])) {
        return (vec3){0.0/0.0, 0.0/0.0, 0.0/0.0};
    }
    return mat3x3_mul_vec3(inv, v);
}

#endif
