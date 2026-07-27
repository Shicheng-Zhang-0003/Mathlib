#include "ml_compiler.h"
#include "ml_optimization.h"

/* v11S CLOSURE IP-15: optimization robustness */

ML_API double ml_optimize_golden(ml_opt_func_t f, double a, double b, double tol, int max_iter) {
    if (ML_UNLIKELY(f == NULL)) {
        return ml_make_nan();
    }

    if (ML_UNLIKELY(tol <= 0.0 || max_iter <= 0)) {
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

    double phi = (1.0 + ml_sqrt(5.0)) / 2.0;
    double resphi = 2.0 - phi;

    double x1 = a + resphi * (b - a);
    double x2 = b - resphi * (b - a);

    double f1 = f(x1);
    double f2 = f(x2);

    if (ML_UNLIKELY(ml_isnan(f1) || ml_isnan(f2))) {
        return ml_make_nan();
    }

    if (ML_UNLIKELY(!ml_isfinite(f1) || !ml_isfinite(f2))) {
        return ml_make_nan();
    }

    for (int i = 0; i < max_iter; i++) {
        if (b - a < tol) {
            break;
        }

        if (ML_UNLIKELY(ml_isnan(f1) || ml_isnan(f2))) {
            return ml_make_nan();
        }

        if (f1 < f2) {
            b = x2;
            x2 = x1;
            f2 = f1;

            x1 = a + resphi * (b - a);
            f1 = f(x1);
        } else {
            a = x1;
            x1 = x2;
            f1 = f2;

            x2 = b - resphi * (b - a);
            f2 = f(x2);
        }

        if (ML_UNLIKELY(x1 == x2)) {
            break;
        }
    }

    double result = (a + b) * 0.5;

    if (ML_UNLIKELY(ml_isnan(result))) {
        return ml_make_nan();
    }

    return result;
}

ML_API double ml_optimize_gradient_descent(ml_opt_func_t f, double start, double lr, double tol, int max_iter) {
    if (ML_UNLIKELY(f == NULL)) {
        return ml_make_nan();
    }

    if (ML_UNLIKELY(lr <= 0.0 || tol <= 0.0 || max_iter <= 0)) {
        return ml_make_nan();
    }

    if (ML_UNLIKELY(ml_isnan(start))) {
        return ml_make_nan();
    }

    double x = start;

    for (int i = 0; i < max_iter; i++) {
        if (ML_UNLIKELY(!ml_isfinite(x))) {
            return ml_make_nan();
        }

        double grad = ml_derivative(f, x, 1e-5);

        if (ML_UNLIKELY(ml_isnan(grad))) {
            return ml_make_nan();
        }

        double x_new = x - lr * grad;

        if (ML_UNLIKELY(ml_isnan(x_new))) {
            return ml_make_nan();
        }

        if (ml_fabs(x_new - x) < tol) {
            return x_new;
        }

        x = x_new;
    }

    return x;
}
