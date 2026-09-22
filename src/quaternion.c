#include "ml_compiler.h"
#include "ml_quaternion.h"

/* v11S CLOSURE IP-16: quaternion edge hardening */

static ml_quat ml_quat_nan(void) {
    double n = ml_make_nan();
    return (ml_quat){n, n, n, n};
}

static int ml_quat_is_finite(ml_quat q) {
    return ml_isfinite(q.w) &&
           ml_isfinite(q.x) &&
           ml_isfinite(q.y) &&
           ml_isfinite(q.z);
}

static ml_quat ml_quat_lerp(ml_quat a, ml_quat b, double t) {
    return (ml_quat){
        a.w + t * (b.w - a.w),
        a.x + t * (b.x - a.x),
        a.y + t * (b.y - a.y),
        a.z + t * (b.z - a.z)
    };
}

ML_API ml_quat ml_quat_mul(ml_quat a, ml_quat b) {
    return (ml_quat){
        a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z,
        a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y,
        a.w*b.y - a.x*b.z + a.y*b.w + a.z*b.x,
        a.w*b.z + a.x*b.y - a.y*b.x + a.z*b.w
    };
}

ML_API ml_quat ml_quat_conjugate(ml_quat q) {
    return (ml_quat){q.w, -q.x, -q.y, -q.z};
}

ML_API ml_quat ml_quat_inverse(ml_quat q) {
    if (!ml_quat_is_finite(q)) return ml_quat_nan();
    {
        /* Max-scaled norm (same as normalize): unscaled w^2+x^2+y^2+z^2
         * overflows for 1e200 components though the inverse exists. */
        double aw = ml_fabs(q.w), ax = ml_fabs(q.x);
        double ay = ml_fabs(q.y), az = ml_fabs(q.z);
        double mx = aw;
        if (ax > mx) mx = ax;
        if (ay > mx) mx = ay;
        if (az > mx) mx = az;
        if (mx == 0.0) return ml_quat_nan();
        {
            double sw = q.w / mx, sx = q.x / mx;
            double sy = q.y / mx, sz = q.z / mx;
            double n2 = sw*sw + sx*sx + sy*sy + sz*sz;
            if (n2 == 0.0 || !ml_isfinite(n2)) return ml_quat_nan();
            /* inv = conj/(mx^2*n2): fold one mx into each component. */
            double inv = 1.0 / (mx * n2);
            return (ml_quat){sw*inv, -sx*inv, -sy*inv, -sz*inv};
        }
    }
}

ML_API ml_quat ml_quat_from_axis_angle(double ax, double ay, double az, double angle) {
    if (!ml_isfinite(ax) || !ml_isfinite(ay) || !ml_isfinite(az) || !ml_isfinite(angle)) {
        return ml_quat_nan();
    }
    {
        /* Scaled axis norm: 1e200 axis overflowed sqrt, 1e-200
         * underflowed to 0->NaN though direction is well-defined. */
        double fx = ml_fabs(ax), fy = ml_fabs(ay), fz = ml_fabs(az);
        double mx = fx;
        if (fy > mx) mx = fy;
        if (fz > mx) mx = fz;
        if (mx == 0.0 || !ml_isfinite(mx)) return ml_quat_nan();
        {
            double sx = ax / mx, sy = ay / mx, sz = az / mx;
            double n = mx * ml_sqrt(sx*sx + sy*sy + sz*sz);
            if (n == 0.0 || !ml_isfinite(n)) return ml_quat_nan();
            {
                double h = angle * 0.5;
                double s = ml_sin(h) / n;
                return ml_quat_normalize((ml_quat){ml_cos(h), ax*s, ay*s, az*s});
            }
        }
    }
}

/*
 * Robust quaternion normalization.
 *
 * We scale by the largest component before computing the norm.
 * This prevents:
 *   - overflow for very large quaternions
 *   - underflow for very small quaternions
 *
 * A zero quaternion has no valid normalized direction, so we return NaN.
 */
ML_API ml_quat ml_quat_normalize(ml_quat q) {
    if (!ml_quat_is_finite(q)) {
        return ml_quat_nan();
    }

    double aw = ml_fabs(q.w);
    double ax = ml_fabs(q.x);
    double ay = ml_fabs(q.y);
    double az = ml_fabs(q.z);

    double max = aw;
    if (ax > max) max = ax;
    if (ay > max) max = ay;
    if (az > max) max = az;

    if (max == 0.0) {
        return ml_quat_nan();
    }

    double sw = q.w / max;
    double sx = q.x / max;
    double sy = q.y / max;
    double sz = q.z / max;

    double mag2 = sw*sw + sx*sx + sy*sy + sz*sz;

    if (mag2 == 0.0 || ml_isnan(mag2)) {
        return ml_quat_nan();
    }

    double inv = 1.0 / ml_sqrt(mag2);

    return (ml_quat){
        sw * inv,
        sx * inv,
        sy * inv,
        sz * inv
    };
}

/*
 * Robust spherical linear interpolation.
 *
 * Contract:
 *   - NaN / Inf inputs produce NaN quaternion
 *   - non-unit inputs are normalized first
 *   - shortest-path interpolation is used (dot < 0 flips b)
 *   - near-parallel quaternions fall back to normalized lerp
 *   - division by sin(theta) is guarded
 */
ML_API ml_quat ml_quat_slerp(ml_quat a, ml_quat b, double t) {
    if (!ml_isfinite(t)) {
        return ml_quat_nan();
    }

    ml_quat an = ml_quat_normalize(a);
    ml_quat bn = ml_quat_normalize(b);

    if (!ml_quat_is_finite(an) || !ml_quat_is_finite(bn)) {
        return ml_quat_nan();
    }

    if (t == 0.0) {
        return an;
    }

    if (t == 1.0) {
        return bn;
    }

    double dot = an.w*bn.w + an.x*bn.x + an.y*bn.y + an.z*bn.z;

    if (!ml_isfinite(dot)) {
        return ml_quat_nan();
    }

    /*
     * Shortest-path correction.
     * If dot is negative, interpolate toward -b instead of b.
     */
    if (dot < 0.0) {
        bn.w = -bn.w;
        bn.x = -bn.x;
        bn.y = -bn.y;
        bn.z = -bn.z;
        dot = -dot;
    }

    /*
     * Clamp into valid acos domain.
     * Rounding can push dot slightly outside [-1, 1]. Note dot>=0 here
     * because the shortest-path flip above already negated negative dots,
     * so the lower clamp is defensive only.
     */
    if (dot > 1.0) dot = 1.0;
    if (dot < 0.0) dot = 0.0;

    /*
     * Near-parallel fallback.
     *
     * When dot is extremely close to 1, sin(theta) becomes numerically
     * unstable. Normalized lerp is the standard safe fallback.
     */
    if (dot > 0.9995) {
        return ml_quat_normalize(ml_quat_lerp(an, bn, t));
    }

    double theta = ml_acos(dot);
    double sin_theta = ml_sin(theta);

    if (!ml_isfinite(sin_theta) || ml_fabs(sin_theta) < 1e-12) {
        return ml_quat_normalize(ml_quat_lerp(an, bn, t));
    }

    double s0 = ml_sin((1.0 - t) * theta) / sin_theta;
    double s1 = ml_sin(t * theta) / sin_theta;

    if (!ml_isfinite(s0) || !ml_isfinite(s1)) {
        return ml_quat_normalize(ml_quat_lerp(an, bn, t));
    }

    /* Renormalize: s0*a+s1*b drifts ~1e-16 per slerp; without this,
     * chained interpolations accumulate norm error. */
    return ml_quat_normalize((ml_quat){
        s0 * an.w + s1 * bn.w,
        s0 * an.x + s1 * bn.x,
        s0 * an.y + s1 * bn.y,
        s0 * an.z + s1 * bn.z
    });
}
