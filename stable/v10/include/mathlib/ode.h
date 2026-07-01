#ifndef LIBMATHC_ODE_H
#define LIBMATHC_ODE_H

typedef double (*ode_func)(double t, double y);

static inline double ode_euler(ode_func f, double t0, double y0, double dt, int steps) {
    double t = t0, y = y0;
    for (int i = 0; i < steps; i++) {
        y += dt * f(t, y);
        t += dt;
    }
    return y;
}

static inline double ode_rk4(ode_func f, double t0, double y0, double dt, int steps) {
    double t = t0, y = y0;
    for (int i = 0; i < steps; i++) {
        double k1 = f(t, y);
        double k2 = f(t + 0.5 * dt, y + 0.5 * dt * k1);
        double k3 = f(t + 0.5 * dt, y + 0.5 * dt * k2);
        double k4 = f(t + dt, y + dt * k3);
        y += (dt / 6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
        t += dt;
    }
    return y;
}
#endif
