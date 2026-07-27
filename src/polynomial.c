#include "ml_compiler.h"
#include "ml_polynomial.h"

/* v11S CLOSURE IP-15: polynomial robustness */

ML_API double ml_polynomial_eval(const double *coeffs, int degree, double x) {
    if (ML_UNLIKELY(coeffs == NULL || degree < 0)) {
        return ml_make_nan();
    }

    double result = coeffs[degree];

    for (int i = degree - 1; i >= 0; i--) {
        result = result * x + coeffs[i];
    }

    return result;
}

ML_API void ml_polynomial_derivative(const double *coeffs, int degree, double *out) {
    if (ML_UNLIKELY(coeffs == NULL || out == NULL || degree < 0)) {
        return;
    }

    for (int i = 0; i < degree; i++) {
        out[i] = coeffs[i + 1] * (double)(i + 1);
    }
}

ML_API double ml_polynomial_newton(const double *coeffs, int degree, double x0, double epsilon, int max_iter) {
    if (ML_UNLIKELY(coeffs == NULL || degree < 0)) {
        return ml_make_nan();
    }

    if (ML_UNLIKELY(epsilon <= 0.0 || max_iter <= 0)) {
        return ml_make_nan();
    }

    if (ML_UNLIKELY(degree == 0)) {
        return ml_make_nan();
    }

    double x = x0;

    for (int iter = 0; iter < max_iter; iter++) {
        double fx = coeffs[degree];
        double dfx = 0.0;

        for (int i = degree - 1; i >= 0; i--) {
            dfx = dfx * x + fx;
            fx = fx * x + coeffs[i];
        }

        if (ML_UNLIKELY(ml_isnan(fx) || ml_isnan(dfx))) {
            return ml_make_nan();
        }

        if (fx == 0.0) {
            return x;
        }

        if (ML_UNLIKELY(ml_fabs(dfx) < epsilon)) {
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
