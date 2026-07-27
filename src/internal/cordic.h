#ifndef LIBMATHC_CORDIC_H
#define LIBMATHC_CORDIC_H
#include "ml_core.h"

static const double cordic_atan[] = {
    0.7853981633974483, 0.4636476090008061, 0.24497866312686414,
    0.12435499454676144, 0.06241880999595735, 0.031239833430268277,
    0.01562372862047683, 0.00781234106010111, 0.003906230131966971,
    0.001953122516478818, 0.000976562189559319, 0.000488281211194898,
    0.000244140620149362, 0.000122070311893670, 0.000061035156174208,
    0.000030517578115521, 0.000015258789061315, 0.000007629394531101,
    0.000003814697265606, 0.000001907348632810, 0.000000953674316405,
    0.000000476837158203, 0.000000238418579101, 0.000000119209289550
};

#define CORDIC_GAIN 0.607252935008881

static inline void ml_cordic_sincos(double theta, double *sin_out, double *cos_out) {
    /* MATHLIB_CLOSURE_P2_P0_5_CORDIC_NONFINITE_GUARD */
    if (ML_UNLIKELY(sin_out == NULL || cos_out == NULL)) {
        return;
    }

    /*
     * Non-finite inputs must propagate as NaN.
     *
     * The old code allowed NaN to fall through the CORDIC iteration,
     * because comparisons against NaN are false, producing finite-looking
     * garbage instead of NaN.
     */
    if (ml_isnan(theta) || ml_isinf(theta)) {
        *sin_out = ml_make_nan();
        *cos_out = ml_make_nan();
        return;
    }

    /* Range reduction to [-pi, pi] */
    theta = ml_fmod(theta, 2.0 * ML_PI);

    /*
     * Defensive guard:
     * finite input should not produce non-finite reduction, but if it does,
     * fail loudly as NaN instead of casting or iterating on garbage.
     */
    if (ml_isnan(theta) || ml_isinf(theta)) {
        *sin_out = ml_make_nan();
        *cos_out = ml_make_nan();
        return;
    }

    if (theta > ML_PI) theta -= 2.0 * ML_PI;
    if (theta < -ML_PI) theta += 2.0 * ML_PI;

    /* Quadrant mapping: CORDIC only converges in [-pi/2, pi/2] */
    int negate_cos = 0;

    if (theta > ML_PI / 2.0) {
        theta = ML_PI - theta;
        negate_cos = 1;
    } else if (theta < -ML_PI / 2.0) {
        theta = -ML_PI - theta;
        negate_cos = 1;
    }

    double x = CORDIC_GAIN;
    double y = 0.0;
    double z = theta;

    for (int i = 0; i < 24; i++) {
        double x_new, y_new;

        if (z >= 0) {
            x_new = x - (y / (double)(1LL << i));
            y_new = y + (x / (double)(1LL << i));
            z -= cordic_atan[i];
        } else {
            x_new = x + (y / (double)(1LL << i));
            y_new = y - (x / (double)(1LL << i));
            z += cordic_atan[i];
        }

        x = x_new;
        y = y_new;
    }

    if (negate_cos) x = -x;

    *cos_out = x;
    *sin_out = y;
}

#endif
