#include "ml_compiler.h"
#include "ml_quadratics.h"
#include <float.h>

/* v11S CLOSURE IP-15: quadratics robustness */

static int ml_quad_nonfinite(double v) {
    return ml_isnan(v) || ml_isinf(v);
}

/* Cast long-double root to double with overflow saturation.
 * Avoids UB of out-of-range (long double)->double conversion. */
static double ml_quad_cast_ld(long double r) {
    if (!(r == r)) {
        return ml_make_nan();
    }
    if (r > (long double)DBL_MAX) {
        return ml_make_inf(0);
    }
    if (r < -(long double)DBL_MAX) {
        return ml_make_inf(1);
    }
    return (double)r;
}

ML_API double ml_equation(double a, double b, double c, double x) {
    /* (a*x + b)*x + c with FMA: 2 roundings instead of 4. */
    return ML_FMA(ML_FMA(a, x, b), x, c);
}

ML_API double ml_formula_pos(double a, double b, double c) {
    if (ml_quad_nonfinite(a) || ml_quad_nonfinite(b) || ml_quad_nonfinite(c)) {
        return ml_make_nan();
    }

    /*
     * Degenerate linear case:
     * bx + c = 0
     */
    if (a == 0.0) {
        if (b == 0.0) {
            return ml_make_nan();
        }

        return -c / b;
    }

    /* Discriminant ALWAYS in long double: double b*b-4*a*c can
     * overflow to Inf or round a tiny positive to 0/negative even
     * though the true value is finite and positive. */
    {
        long double disc_ld = (long double)b * (long double)b - 4.0L * (long double)a * (long double)c;
        long double sqrt_ld = 0.0L;
        double sqrt_disc = 0.0;
        if (!(disc_ld >= 0.0L)) {
            return ml_make_nan();
        }
        if (disc_ld == 0.0L) {
            return -b / (2.0 * a);
        }
        sqrt_ld = __builtin_sqrtl(disc_ld);
        sqrt_disc = (double)sqrt_ld;
        if (!ml_isfinite(sqrt_disc) || (sqrt_disc == 0.0 && sqrt_ld != 0.0L)) {
            /* Cast underflowed to 0 or overflowed to Inf: compute the
             * stable q-form in long double, then saturate on cast. */
            long double bl = (long double)b;
            long double al = (long double)a;
            long double cl = (long double)c;
            long double q_ld = 0.0L;
            if (bl >= 0.0L) {
                q_ld = -0.5L * (bl + sqrt_ld);
                if (q_ld == 0.0L) {
                    return -b / (2.0 * a);
                }
                return ml_quad_cast_ld(cl / q_ld);
            }
            q_ld = -0.5L * (bl - sqrt_ld);
            return ml_quad_cast_ld(q_ld / al);
        }
        /*
         * Numerically stable branch:
         * Use q-form to avoid catastrophic cancellation.
         */
        if (b >= 0.0) {
            double q = -0.5 * (b + sqrt_disc);

            if (q == 0.0) {
                return -b / (2.0 * a);
            }

            return c / q;
        } else {
            double q = -0.5 * (b - sqrt_disc);
            return q / a;
        }
    }
}

ML_API double ml_formula_neg(double a, double b, double c) {
    if (ml_quad_nonfinite(a) || ml_quad_nonfinite(b) || ml_quad_nonfinite(c)) {
        return ml_make_nan();
    }

    /*
     * Degenerate linear case:
     * bx + c = 0
     */
    if (a == 0.0) {
        if (b == 0.0) {
            return ml_make_nan();
        }

        return -c / b;
    }

    /* Discriminant ALWAYS in long double (see ml_formula_pos). */
    {
        long double disc_ld = (long double)b * (long double)b - 4.0L * (long double)a * (long double)c;
        long double sqrt_ld = 0.0L;
        double sqrt_disc = 0.0;
        if (!(disc_ld >= 0.0L)) {
            return ml_make_nan();
        }
        if (disc_ld == 0.0L) {
            return -b / (2.0 * a);
        }
        sqrt_ld = __builtin_sqrtl(disc_ld);
        sqrt_disc = (double)sqrt_ld;
        if (!ml_isfinite(sqrt_disc) || (sqrt_disc == 0.0 && sqrt_ld != 0.0L)) {
            /* Cast underflow/overflow: long-double q-form + saturation. */
            long double bl = (long double)b;
            long double al = (long double)a;
            long double cl = (long double)c;
            long double q_ld = 0.0L;
            if (bl >= 0.0L) {
                q_ld = -0.5L * (bl + sqrt_ld);
                return ml_quad_cast_ld(q_ld / al);
            }
            q_ld = -0.5L * (bl - sqrt_ld);
            if (q_ld == 0.0L) {
                return -b / (2.0 * a);
            }
            return ml_quad_cast_ld(cl / q_ld);
        }
        if (b >= 0.0) {
            double q = -0.5 * (b + sqrt_disc);
            return q / a;
        } else {
            double q = -0.5 * (b - sqrt_disc);

            if (q == 0.0) {
                return -b / (2.0 * a);
            }

            return c / q;
        }
    }
}

ML_API int ml_cubic(double a, double b, double c, double d, double *r0, double *r1, double *r2) {
    /* Depressed cubic t^3+pt+q=0 via Cardano with trig fallback for
     * casus irreducibilis (3 real roots). Returns #real roots (0/1/3).
     * a==0 delegates to quadratic. All roots real-valued (complex pr. NaN). */
    if (ml_quad_nonfinite(a) || ml_quad_nonfinite(b) || ml_quad_nonfinite(c) || ml_quad_nonfinite(d)) {
        return 0;
    }
    if (a == 0.0) {
        if (b == 0.0) {
            if (c == 0.0) return 0;
            if (r0) *r0 = -d / c;
            return 1;
        }
        if (r0 && r1) {
            *r0 = ml_formula_pos(b, c, d);
            *r1 = ml_formula_neg(b, c, d);
            if (ml_isnan(*r0)) return 0;
            return (*r0 == *r1) ? 1 : 2;
        }
        return 0;
    }
    if (!r0 || !r1 || !r2) return 0;
    {
        long double inv_a = 1.0L / (long double)a;
        long double B = (long double)b * inv_a / 3.0L;
        long double C = (long double)c * inv_a;
        long double D = (long double)d * inv_a;
        long double p = C - 3.0L * B * B;
        long double q = 2.0L * B * B * B - B * C + D;
        long double disc = (q * q) / 4.0L + (p * p * p) / 27.0L;
        if (disc > 0.0L) {
            long double sq = __builtin_sqrtl(disc);
            long double u = __builtin_cbrtl(-q / 2.0L + sq);
            long double v = __builtin_cbrtl(-q / 2.0L - sq);
            *r0 = (double)(u + v - B);
            *r1 = ml_make_nan();
            *r2 = ml_make_nan();
            return 1;
        }
        if (disc == 0.0L) {
            long double u = __builtin_cbrtl(-q / 2.0L);
            *r0 = (double)(2.0L * u - B);
            *r1 = (double)(-u - B);
            *r2 = ml_make_nan();
            return (*r0 == *r1) ? 1 : 2;
        }
        {
            long double r = __builtin_sqrtl(-p / 3.0L);
            long double denom = 2.0L * r * r * r;
            long double arg = (denom == 0.0L) ? 0.0L : -q / denom;
            /* Clamp: rounding pushes arg outside [-1,1] near multiple
             * roots -> acosl(NaN) -> NaN roots. */
            if (arg > 1.0L) arg = 1.0L;
            if (arg < -1.0L) arg = -1.0L;
            long double th = __builtin_acosl(arg) / 3.0L;
            long double t0 = 2.0L * r * __builtin_cosl(th);
            long double t1 = 2.0L * r * __builtin_cosl(th + 2.0L * 3.14159265358979323846264338327950288L / 3.0L);
            long double t2 = 2.0L * r * __builtin_cosl(th + 4.0L * 3.14159265358979323846264338327950288L / 3.0L);
            *r0 = (double)(t0 - B);
            *r1 = (double)(t1 - B);
            *r2 = (double)(t2 - B);
            return 3;
        }
    }
}
