#include "ml_compiler.h"
#include "ml_ode.h"
#include "ml_core.h"
#include <stddef.h>

/* v11S CLOSURE IP-15: ODE robustness */

ML_API double ml_ode_euler(ml_ode_func_t f, double t0, double y0, double dt, int steps) {
    if (ML_UNLIKELY(f == NULL)) {
        return ml_make_nan();
    }

    if (ML_UNLIKELY(steps < 0)) {
        return ml_make_nan();
    }

    if (ML_UNLIKELY(ml_isnan(dt) || ml_isinf(dt))) {
        return ml_make_nan();
    }

    if (ML_UNLIKELY(ml_isnan(t0) || ml_isnan(y0))) {
        return ml_make_nan();
    }

    if (steps == 0) {
        return y0;
    }

    double t = t0;
    double y = y0;

    for (int i = 0; i < steps; i++) {
        if (ML_UNLIKELY(ml_isnan(t) || ml_isnan(y))) {
            return ml_make_nan();
        }

        double k = f(t, y);

        if (ML_UNLIKELY(ml_isnan(k))) {
            return ml_make_nan();
        }

        y += dt * k;
        t += dt;
    }

    return y;
}

ML_API double ml_ode_rk4(ml_ode_func_t f, double t0, double y0, double dt, int steps) {
    if (ML_UNLIKELY(f == NULL)) {
        return ml_make_nan();
    }

    if (ML_UNLIKELY(steps < 0)) {
        return ml_make_nan();
    }

    if (ML_UNLIKELY(ml_isnan(dt) || ml_isinf(dt))) {
        return ml_make_nan();
    }

    if (ML_UNLIKELY(ml_isnan(t0) || ml_isnan(y0))) {
        return ml_make_nan();
    }

    if (steps == 0) {
        return y0;
    }

    double t = t0;
    double y = y0;

    for (int i = 0; i < steps; i++) {
        if (ML_UNLIKELY(ml_isnan(t) || ml_isnan(y))) {
            return ml_make_nan();
        }

        double k1 = f(t, y);

        if (ML_UNLIKELY(ml_isnan(k1))) {
            return ml_make_nan();
        }

        double k2 = f(t + 0.5 * dt, y + 0.5 * dt * k1);

        if (ML_UNLIKELY(ml_isnan(k2))) {
            return ml_make_nan();
        }

        double k3 = f(t + 0.5 * dt, y + 0.5 * dt * k2);

        if (ML_UNLIKELY(ml_isnan(k3))) {
            return ml_make_nan();
        }

        double k4 = f(t + dt, y + dt * k3);

        if (ML_UNLIKELY(ml_isnan(k4))) {
            return ml_make_nan();
        }

        y += (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        t += dt;
    }

    return y;
}
