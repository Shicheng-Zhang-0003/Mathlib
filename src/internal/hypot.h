#ifndef ML_INTERNAL_HYPOT_H
#define ML_INTERNAL_HYPOT_H

#include "ml_core.h"

/*
 * Overflow-safe hypotenuse.
 *
 * Computes sqrt(x*x + y*y) without unnecessary overflow/underflow.
 */
static inline double ml_hypot_internal(double x, double y) {
    if (ml_isinf(x) || ml_isinf(y)) return ml_make_inf(0);
    if (ml_isnan(x) || ml_isnan(y)) return ml_make_nan();

    double ax = ml_fabs(x);
    double ay = ml_fabs(y);

    if (ax < ay) {
        double t = ax;
        ax = ay;
        ay = t;
    }

    if (ax == 0.0) return 0.0;

    double r = ay / ax;
    return ax * ml_sqrt(1.0 + r * r);
}

#endif
