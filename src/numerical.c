#include "ml_compiler.h"
#include "ml_numerical.h"

/* v11S CLOSURE IP-15: numerical methods robustness */

ML_API double ml_newton_raphson(ml_func_t f, ml_func_t df, double x0, double epsilon, int max_iter) {
    if (ML_UNLIKELY(f == NULL || df == NULL)) {
        return ml_make_nan();
    }

    if (ML_UNLIKELY(epsilon <= 0.0 || max_iter <= 0)) {
        return ml_make_nan();
    }

    if (ML_UNLIKELY(ml_isnan(x0))) {
        return ml_make_nan();
    }

    double x = x0;

    for (int i = 0; i < max_iter; i++) {
        double fx = f(x);

        if (ML_UNLIKELY(ml_isnan(fx))) {
            return ml_make_nan();
        }

        if (fx == 0.0) {
            return x;
        }

        double dfx = df(x);

        if (ML_UNLIKELY(ml_isnan(dfx))) {
            return ml_make_nan();
        }

        if (ML_UNLIKELY(dfx == 0.0 || ml_fabs(dfx) < epsilon)) {
            return ml_make_nan();
        }

        double x_next = x - fx / dfx;

        if (ML_UNLIKELY(ml_isnan(x_next))) {
            return ml_make_nan();
        }

        if (ml_fabs(x_next - x) < epsilon) {
            return x_next;
        }

        x = x_next;
    }

    return x;
}

ML_API double ml_bisection(ml_func_t f, double a, double b, double epsilon, int max_iter) {
    if (ML_UNLIKELY(f == NULL)) {
        return ml_make_nan();
    }

    if (ML_UNLIKELY(epsilon <= 0.0 || max_iter <= 0)) {
        return ml_make_nan();
    }

    if (ML_UNLIKELY(ml_isnan(a) || ml_isnan(b) || ml_isinf(a) || ml_isinf(b))) {
        return ml_make_nan();
    }

    if (a > b) {
        double t = a;
        a = b;
        b = t;
    }

    double fa = f(a);
    double fb = f(b);

    if (ML_UNLIKELY(ml_isnan(fa) || ml_isnan(fb))) {
        return ml_make_nan();
    }

    if (ML_UNLIKELY(!ml_isfinite(fa) || !ml_isfinite(fb))) {
        return ml_make_nan();
    }

    if (fa == 0.0) {
        return a;
    }

    if (fb == 0.0) {
        return b;
    }

    if (ml_signbit(fa) == ml_signbit(fb)) {
        return ml_make_nan();
    }

    double c = a;

    for (int i = 0; i < max_iter; i++) {
        c = (a + b) * 0.5;

        /*
         * If the midpoint stops changing, further bisection is impossible
         * in double precision. Return the best current estimate.
         */
        if (c == a || c == b) {
            return c;
        }

        double fc = f(c);

        if (ML_UNLIKELY(ml_isnan(fc) || !ml_isfinite(fc))) {
            return ml_make_nan();
        }

        if (fc == 0.0 || ml_fabs(b - a) < epsilon) {
            return c;
        }

        if (ml_signbit(fa) != ml_signbit(fc)) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    }

    return c;
}

ML_API double ml_derivative(ml_func_t f, double x, double h) {
    if (ML_UNLIKELY(f == NULL)) {
        return ml_make_nan();
    }

    if (ML_UNLIKELY(h == 0.0 || ml_isnan(h) || ml_isinf(h))) {
        return ml_make_nan();
    }

    if (ML_UNLIKELY(ml_isnan(x) || ml_isinf(x))) {
        return ml_make_nan();
    }

    double xp = x + h;
    double xm = x - h;

    if (ML_UNLIKELY(!ml_isfinite(xp) || !ml_isfinite(xm))) {
        return ml_make_nan();
    }

    if (ML_UNLIKELY(xp == xm)) {
        return ml_make_nan();
    }

    double fp = f(xp);
    double fm = f(xm);

    if (ML_UNLIKELY(ml_isnan(fp) || ml_isnan(fm))) {
        return ml_make_nan();
    }

    return (fp - fm) / (2.0 * h);
}

ML_API double ml_second_derivative(ml_func_t f, double x, double h) {
    if (ML_UNLIKELY(f == NULL)) {
        return ml_make_nan();
    }

    if (ML_UNLIKELY(h == 0.0 || ml_isnan(h) || ml_isinf(h))) {
        return ml_make_nan();
    }

    if (ML_UNLIKELY(ml_isnan(x) || ml_isinf(x))) {
        return ml_make_nan();
    }

    double xp = x + h;
    double xm = x - h;

    if (ML_UNLIKELY(!ml_isfinite(xp) || !ml_isfinite(xm))) {
        return ml_make_nan();
    }

    if (ML_UNLIKELY(xp == xm)) {
        return ml_make_nan();
    }

    double fp = f(xp);
    double f0 = f(x);
    double fm = f(xm);

    if (ML_UNLIKELY(ml_isnan(fp) || ml_isnan(f0) || ml_isnan(fm))) {
        return ml_make_nan();
    }

    double denom = h * h;

    if (ML_UNLIKELY(denom == 0.0 || !ml_isfinite(denom))) {
        return ml_make_nan();
    }

    return (fp - 2.0 * f0 + fm) / denom;
}

ML_API double ml_integral_simpson(ml_func_t f, double a, double b, int n) {
    if (ML_UNLIKELY(f == NULL)) {
        return ml_make_nan();
    }

    /*
     * Simpson's rule strictly requires an even number of subintervals.
     * Reject invalid grids explicitly instead of silently changing them.
     */
    if (ML_UNLIKELY(n < 2 || (n % 2) == 1)) {
        return ml_make_nan();
    }

    if (ML_UNLIKELY(ml_isnan(a) || ml_isnan(b) || ml_isinf(a) || ml_isinf(b))) {
        return ml_make_nan();
    }

    if (a == b) {
        return 0.0;
    }

    double h = (b - a) / (double)n;

    if (ML_UNLIKELY(h == 0.0 || ml_isnan(h) || ml_isinf(h))) {
        return ml_make_nan();
    }

    double fa = f(a);
    double fb = f(b);

    if (ML_UNLIKELY(ml_isnan(fa) || ml_isnan(fb))) {
        return ml_make_nan();
    }

    double result = fa + fb;

    for (int i = 1; i < n; i++) {
        double x = a + (double)i * h;
        double fx = f(x);

        if (ML_UNLIKELY(ml_isnan(fx))) {
            return ml_make_nan();
        }

        if ((i % 2) == 0) {
            result += 2.0 * fx;
        } else {
            result += 4.0 * fx;
        }
    }

    return result * h / 3.0;
}
