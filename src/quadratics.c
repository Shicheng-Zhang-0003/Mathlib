#include "ml_compiler.h"
#include "ml_quadratics.h"

/* v11S CLOSURE IP-15: quadratics robustness */

static int ml_quad_nonfinite(double v) {
    return ml_isnan(v) || ml_isinf(v);
}

ML_API double ml_equation(double a, double b, double c, double x) {
    return (a * x * x) + b * x + c;
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

    double disc = b * b - 4.0 * a * c;

    if (ml_isnan(disc) || ml_isinf(disc)) {
        return ml_make_nan();
    }

    if (disc < 0.0) {
        return ml_make_nan();
    }

    if (disc == 0.0) {
        return -b / (2.0 * a);
    }

    double sqrt_disc = ml_sqrt(disc);

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

    double disc = b * b - 4.0 * a * c;

    if (ml_isnan(disc) || ml_isinf(disc)) {
        return ml_make_nan();
    }

    if (disc < 0.0) {
        return ml_make_nan();
    }

    if (disc == 0.0) {
        return -b / (2.0 * a);
    }

    double sqrt_disc = ml_sqrt(disc);

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
