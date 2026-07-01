#ifndef LIBMATHC_CORDIC_H
#define LIBMATHC_CORDIC_H
#include "ml_core.h"
#include "payne_hanek.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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
    // 1. Range reduction to [-pi, pi]
    theta = ml_reduce_payne_hanek(theta);
    if (theta > M_PI) theta -= 2.0 * M_PI;
    if (theta < -M_PI) theta += 2.0 * M_PI;

    // 2. Quadrant Mapping: CORDIC only converges in [-pi/2, pi/2]
    int negate_cos = 0;
    if (theta > M_PI / 2.0) {
        theta = M_PI - theta;
        negate_cos = 1;
    } else if (theta < -M_PI / 2.0) {
        theta = -M_PI - theta;
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

    // 3. Apply Quadrant Sign Correction
    if (negate_cos) x = -x;

    *cos_out = x;
    *sin_out = y;
}
#endif
