#include "ml_compiler.h"
#include "ml_ode.h"
#include "ml_core.h"
#include "ml_exp_log.h"
#include "ml_trig.h"
#include <stddef.h>

/* v11S CLOSURE IP-15: ODE robustness */

ML_API double ml_ode_euler(ml_ode_func_t f, double t0, double y0, double dt, int steps) {
    if (ML_UNLIKELY(f == NULL)) {
        return ml_make_nan();
    }

    if (ML_UNLIKELY(steps < 0)) {
        return ml_make_nan();
    }

    if (ML_UNLIKELY(!ml_isfinite(dt))) {
        return ml_make_nan();
    }

    if (ML_UNLIKELY(!ml_isfinite(t0) || !ml_isfinite(y0))) {
        return ml_make_nan();
    }

    if (steps == 0) {
        return y0;
    }

    double t = t0;
    double y = y0;

    for (int i = 0; i < steps; i++) {
        if (ML_UNLIKELY(!ml_isfinite(t) || !ml_isfinite(y))) {
            return ml_make_nan();
        }

        double k = f(t, y);

        if (ML_UNLIKELY(!ml_isfinite(k))) {
            return ml_make_nan();
        }

        y += dt * k;
        t += dt;
        if (ML_UNLIKELY(!ml_isfinite(y) || !ml_isfinite(t))) {
            return ml_make_nan();
        }
    }

    if (ML_UNLIKELY(!ml_isfinite(y))) {
        return ml_make_nan();
    }
    return y;
}

ML_API double ml_ode_rk4(ml_ode_func_t f, double t0, double y0, double dt, int steps) {
    if (ML_UNLIKELY(f == NULL)) return ml_make_nan();
    if (ML_UNLIKELY(steps < 0)) return ml_make_nan();
    if (ML_UNLIKELY(!ml_isfinite(dt))) return ml_make_nan();
    if (ML_UNLIKELY(!ml_isfinite(t0) || !ml_isfinite(y0))) return ml_make_nan();
    if (steps == 0) return y0;
    {
        double t = t0, y = y0;
        for (int i = 0; i < steps; i++) {
            double k1, k2, k3, k4;
            if (ML_UNLIKELY(!ml_isfinite(t) || !ml_isfinite(y))) return ml_make_nan();
            k1 = f(t, y);
            if (ML_UNLIKELY(!ml_isfinite(k1))) return ml_make_nan();
            k2 = f(t + 0.5 * dt, y + 0.5 * dt * k1);
            if (ML_UNLIKELY(!ml_isfinite(k2))) return ml_make_nan();
            k3 = f(t + 0.5 * dt, y + 0.5 * dt * k2);
            if (ML_UNLIKELY(!ml_isfinite(k3))) return ml_make_nan();
            k4 = f(t + dt, y + dt * k3);
            if (ML_UNLIKELY(!ml_isfinite(k4))) return ml_make_nan();
            y += (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
            t += dt;
            if (ML_UNLIKELY(!ml_isfinite(y) || !ml_isfinite(t))) return ml_make_nan();
        }
        return y;
    }
}

ML_API double ml_ode_leapfrog(ml_ode_func_t f, double t0, double y0, double dt, int steps) {
    /* Midpoint RK2 (explicit midpoint): y_{n+1} = y_n + dt*f(t+dt/2,
     * y+dt/2*f(t,y)). This is second-order but NOT symplectic
     * velocity-Verlet (which needs a 2nd-order Hamiltonian split).
     * Do not use as an energy-preserving integrator. */
    if (ML_UNLIKELY(f == NULL)) return ml_make_nan();
    if (ML_UNLIKELY(steps < 0)) return ml_make_nan();
    if (ML_UNLIKELY(!ml_isfinite(dt))) return ml_make_nan();
    if (ML_UNLIKELY(!ml_isfinite(t0) || !ml_isfinite(y0))) return ml_make_nan();
    if (steps == 0) return y0;
    {
        double t = t0, y = y0;
        for (int i = 0; i < steps; i++) {
            double k1, k2;
            if (ML_UNLIKELY(!ml_isfinite(t) || !ml_isfinite(y))) return ml_make_nan();
            k1 = f(t, y);
            if (ML_UNLIKELY(!ml_isfinite(k1))) return ml_make_nan();
            {
                double yh = y + 0.5 * dt * k1;
                double th = t + 0.5 * dt;
                k2 = f(th, yh);
                if (ML_UNLIKELY(!ml_isfinite(k2))) return ml_make_nan();
                y += dt * k2;
                t += dt;
            }
            if (ML_UNLIKELY(!ml_isfinite(y) || !ml_isfinite(t))) return ml_make_nan();
        }
        return y;
    }
}

ML_API double ml_ode_dp5(ml_ode_func_t f, double t0, double y0, double t1, double rtol, double atol) {
    /* Dormand-Prince 5(4) adaptive with I-controller. Scalar only.
     * rtol/atol must be finite positive; t1 must be finite. */
    if (ML_UNLIKELY(f == NULL)) return ml_make_nan();
    if (ML_UNLIKELY(!(rtol > 0.0) || !ml_isfinite(rtol))) return ml_make_nan();
    if (ML_UNLIKELY(!(atol > 0.0) || !ml_isfinite(atol))) return ml_make_nan();
    if (ML_UNLIKELY(!ml_isfinite(t0) || !ml_isfinite(y0) || !ml_isfinite(t1))) return ml_make_nan();
    if (t1 == t0) return y0;
    {
        static const double C2 = 1.0/5.0, C3 = 3.0/10.0, C4 = 4.0/5.0,
            C5 = 8.0/9.0, C6 = 1.0, C7 = 1.0;
        static const double A21 = 1.0/5.0;
        static const double A31 = 3.0/40.0, A32 = 9.0/40.0;
        static const double A41 = 44.0/45.0, A42 = -56.0/15.0, A43 = 32.0/9.0;
        static const double A51 = 19372.0/6561.0, A52 = -25360.0/2187.0,
            A53 = 64448.0/6561.0, A54 = -212.0/729.0;
        static const double A61 = 9017.0/3168.0, A62 = -355.0/33.0,
            A63 = 46732.0/5247.0, A64 = 49.0/176.0, A65 = -5103.0/18656.0;
        static const double A71 = 35.0/384.0, A73 = 500.0/1113.0,
            A74 = 125.0/192.0, A75 = -2187.0/6784.0, A76 = 11.0/84.0;
        static const double E1 = 71.0/57600.0, E3 = -71.0/16695.0,
            E4 = 71.0/1920.0, E5 = -17253.0/339200.0,
            E6 = 22.0/525.0, E7 = -1.0/40.0;
        (void)C2; (void)C3; (void)C4; (void)C5; (void)C6; (void)C7;
        double t = t0, y = y0;
        double h = (t1 - t0) * 0.01;
        if (!ml_isfinite(h) || h == 0.0) return ml_make_nan();
        if ((t1 > t0 && h < 0.0) || (t1 < t0 && h > 0.0)) h = -h;
        for (int iter = 0; iter < 100000; iter++) {
            double k1, k2, k3, k4, k5, k6, k7, y5, err;
            if ((h > 0.0 && t + h > t1) || (h < 0.0 && t + h < t1)) h = t1 - t;
            k1 = f(t, y);
            if (!ml_isfinite(k1)) return ml_make_nan();
            k2 = f(t + h*(1.0/5.0), y + h*(A21*k1));
            if (!ml_isfinite(k2)) return ml_make_nan();
            k3 = f(t + h*(3.0/10.0), y + h*(A31*k1 + A32*k2));
            if (!ml_isfinite(k3)) return ml_make_nan();
            k4 = f(t + h*(4.0/5.0), y + h*(A41*k1 + A42*k2 + A43*k3));
            if (!ml_isfinite(k4)) return ml_make_nan();
            k5 = f(t + h*(8.0/9.0), y + h*(A51*k1 + A52*k2 + A53*k3 + A54*k4));
            if (!ml_isfinite(k5)) return ml_make_nan();
            k6 = f(t + h, y + h*(A61*k1 + A62*k2 + A63*k3 + A64*k4 + A65*k5));
            if (!ml_isfinite(k6)) return ml_make_nan();
            y5 = y + h*(A71*k1 + A73*k3 + A74*k4 + A75*k5 + A76*k6);
            k7 = f(t + h, y5);
            if (!ml_isfinite(k7)) return ml_make_nan();
            err = h*(E1*k1 + E3*k3 + E4*k4 + E5*k5 + E6*k6 + E7*k7);
            {
                double sc = atol + rtol * (ml_fabs(y) > ml_fabs(y5) ? ml_fabs(y) : ml_fabs(y5));
                double enorm = ml_fabs(err) / sc;
                if (enorm <= 1.0) {
                    t += h; y = y5;
                    if (t == t1) return y;
                }
                {
                    double fac = (enorm == 0.0) ? 5.0 : 0.9 / ml_pow(enorm, 0.2);
                    if (fac < 0.2) fac = 0.2;
                    if (fac > 5.0) fac = 5.0;
                    h *= fac;
                    if (!ml_isfinite(h) || h == 0.0) return ml_make_nan();
                }
            }
            if (t == t1) return y;
        }
        return ml_make_nan();
    }
}

ML_API double ml_ode_heun(ml_ode_func_t f, double t0, double y0, double dt, int steps) {
    if (ML_UNLIKELY(f == NULL)) return ml_make_nan();
    if (ML_UNLIKELY(steps < 0)) return ml_make_nan();
    if (ML_UNLIKELY(!ml_isfinite(dt))) return ml_make_nan();
    if (ML_UNLIKELY(!ml_isfinite(t0) || !ml_isfinite(y0))) return ml_make_nan();
    if (steps == 0) return y0;
    {
        double t = t0, y = y0;
        for (int i = 0; i < steps; i++) {
            double k1, kp, k2;
            if (ML_UNLIKELY(!ml_isfinite(t) || !ml_isfinite(y))) return ml_make_nan();
            k1 = f(t, y);
            if (ML_UNLIKELY(!ml_isfinite(k1))) return ml_make_nan();
            kp = y + dt * k1;
            k2 = f(t + dt, kp);
            if (ML_UNLIKELY(!ml_isfinite(k2))) return ml_make_nan();
            y += dt * 0.5 * (k1 + k2);
            t += dt;
            if (ML_UNLIKELY(!ml_isfinite(y) || !ml_isfinite(t))) return ml_make_nan();
        }
        return y;
    }
}

ML_API double ml_suvat_s(double u, double a, double t) {
    if (ML_UNLIKELY(!ml_isfinite(u) || !ml_isfinite(a) || !ml_isfinite(t))) return ml_make_nan();
    return u * t + 0.5 * a * t * t;
}

ML_API double ml_suvat_v(double u, double a, double t) {
    if (ML_UNLIKELY(!ml_isfinite(u) || !ml_isfinite(a) || !ml_isfinite(t))) return ml_make_nan();
    return u + a * t;
}

ML_API double ml_projectile_range(double v0, double theta, double g) {
    if (ML_UNLIKELY(!ml_isfinite(v0) || !ml_isfinite(theta) || !(g > 0.0) || !ml_isfinite(g))) return ml_make_nan();
    if (ML_UNLIKELY(v0 < 0.0)) return ml_make_nan();
    {
        double s = ml_sin(2.0 * theta);
        /* Negative sin2θ (θ outside [0,pi/2]) is unphysical range. */
        if (s < 0.0) return ml_make_nan();
        return v0 * v0 * s / g;
    }
}

ML_API double ml_projectile_height(double v0, double theta, double g) {
    if (ML_UNLIKELY(!ml_isfinite(v0) || !ml_isfinite(theta) || !(g > 0.0) || !ml_isfinite(g))) return ml_make_nan();
    if (ML_UNLIKELY(v0 < 0.0)) return ml_make_nan();
    {
        double s = ml_sin(theta);
        return v0 * v0 * s * s / (2.0 * g);
    }
}

ML_API double ml_shm_x(double A, double omega, double t, double phi) {
    if (ML_UNLIKELY(!ml_isfinite(A) || !ml_isfinite(omega) || !ml_isfinite(t) || !ml_isfinite(phi))) {
        return ml_make_nan();
    }
    return A * ml_cos(omega * t + phi);
}

ML_API double ml_collision_1d(double m1, double m2, double u1, double u2, double *v1, double *v2) {
    if (ML_UNLIKELY(!v1 || !v2)) return ml_make_nan();
    if (ML_UNLIKELY(!(m1 > 0.0) || !(m2 > 0.0))) return ml_make_nan();
    if (ML_UNLIKELY(!ml_isfinite(m1) || !ml_isfinite(m2) || !ml_isfinite(u1) || !ml_isfinite(u2))) {
        *v1 = ml_make_nan(); *v2 = ml_make_nan(); return ml_make_nan();
    }
    {
        double M = m1 + m2;
        if (ML_UNLIKELY(!ml_isfinite(M) || M == 0.0)) {
            *v1 = ml_make_nan(); *v2 = ml_make_nan(); return ml_make_nan();
        }
        *v1 = ((m1 - m2) * u1 + 2.0 * m2 * u2) / M;
        *v2 = ((m2 - m1) * u2 + 2.0 * m1 * u1) / M;
        if (ML_UNLIKELY(!ml_isfinite(*v1) || !ml_isfinite(*v2))) {
            *v1 = ml_make_nan(); *v2 = ml_make_nan(); return ml_make_nan();
        }
        return 0.0;
    }
}

ML_API int ml_ode2_const(double a, double b, double c, double *r0, double *r1) {
    if (ML_UNLIKELY(!r0 || !r1)) return 0;
    if (ML_UNLIKELY(!ml_isfinite(a) || !ml_isfinite(b) || !ml_isfinite(c))) return 0;
    if (a == 0.0) {
        if (b == 0.0) return 0;
        *r0 = -c / b; *r1 = ml_make_nan();
        return 1;
    }
    {
        /* Long-double discriminant: naive b*b-4ac overflows to Inf for
         * 1e200-scale coefficients though roots are finite. */
        long double disc = (long double)b * (long double)b - 4.0L * (long double)a * (long double)c;
        if (!(disc >= 0.0L) || !(disc < (long double)1e4932L)) return 0;
        if (disc == 0.0L) { *r0 = -b / (2.0 * a); *r1 = ml_make_nan(); return 1; }
        {
            double s = (double)__builtin_sqrtl(disc);
            if (!ml_isfinite(s)) return 0;
            /* Stable pairing: avoid (-b+s) cancellation. */
            if (b >= 0.0) {
                double q = -0.5 * (b + s);
                if (q == 0.0) { *r0 = -b / (2.0 * a); *r1 = ml_make_nan(); return 1; }
                *r0 = c / q; *r1 = q / a;
            } else {
                double q = -0.5 * (b - s);
                *r0 = q / a; *r1 = c / q;
            }
            return 2;
        }
    }
}
