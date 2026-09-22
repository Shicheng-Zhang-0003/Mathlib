#include "ml_compiler.h"
#include "ml_transforms.h"

ML_API void ml_dct2(const double *x, double *X, int n) {
    if (ML_UNLIKELY(!x || !X || n <= 0)) return;
    for (int k = 0; k < n; k++) {
        long double s = 0.0L;
        for (int i = 0; i < n; i++) {
            if (ML_UNLIKELY(!ml_isfinite(x[i]))) { X[k] = ml_make_nan(); break; }
            long double ang = 3.14159265358979323846264338327950288L
                * (long double)k * (2.0L * (long double)i + 1.0L)
                / (2.0L * (long double)n);
            s += (long double)x[i] * __builtin_cosl(ang);
        }
        X[k] = (double)s;
    }
}

ML_API void ml_dct3(const double *X, double *x, int n) {
    if (ML_UNLIKELY(!X || !x || n <= 0)) return;
    for (int i = 0; i < n; i++) {
        long double s = 0.5L * (long double)X[0];
        for (int k = 1; k < n; k++) {
            if (ML_UNLIKELY(!ml_isfinite(X[k]))) { x[i] = ml_make_nan(); break; }
            long double ang = 3.14159265358979323846264338327950288L
                * (long double)k * (2.0L * (long double)i + 1.0L)
                / (2.0L * (long double)n);
            s += (long double)X[k] * __builtin_cosl(ang);
        }
        x[i] = (double)(s * (2.0L / (long double)n));
    }
}

ML_API void ml_dst2(const double *x, double *X, int n) {
    if (ML_UNLIKELY(!x || !X || n <= 0)) return;
    for (int k = 1; k <= n; k++) {
        long double s = 0.0L;
        for (int i = 1; i <= n; i++) {
            if (ML_UNLIKELY(!ml_isfinite(x[i-1]))) { X[k-1] = ml_make_nan(); break; }
            long double ang = 3.14159265358979323846264338327950288L
                * (long double)k * (long double)i / ((long double)n + 1.0L);
            s += (long double)x[i-1] * __builtin_sinl(ang);
        }
        X[k-1] = (double)s;
    }
}

ML_API double ml_parseval_energy(const double *x, int n) {
    if (ML_UNLIKELY(!x || n <= 0)) return ml_make_nan();
    {
        long double s = 0.0L;
        for (int i = 0; i < n; i++) {
            if (ML_UNLIKELY(!ml_isfinite(x[i]))) return ml_make_nan();
            s += (long double)x[i] * (long double)x[i];
        }
        double r = (double)s;
        return ml_isfinite(r) ? r : ml_make_nan();
    }
}

ML_API ml_status_t ml_circular_conv(const double *a, const double *b, double *out, int n) {
    if (ML_UNLIKELY(!a || !b || !out || n <= 0 || n > 4096)) return ML_ERR_INVALID_ARG;
    for (int i = 0; i < n; i++) {
        if (ML_UNLIKELY(!ml_isfinite(a[i]) || !ml_isfinite(b[i]))) return ML_ERR_NAN_INPUT;
        long double s = 0.0L;
        for (int j = 0; j < n; j++) {
            int k = i - j;
            if (k < 0) k += n;
            s += (long double)a[j] * (long double)b[k];
        }
        out[i] = (double)s;
        if (ML_UNLIKELY(!ml_isfinite(out[i]))) return ML_ERR_INTERNAL;
    }
    return ML_SUCCESS;
}

ML_API double ml_laplace_exp(double s, double a) {
    /* L{e^{at}} = 1/(s-a), s>a. Undergrad Laplace table. */
    if (ML_UNLIKELY(ml_isnan(s) || ml_isnan(a))) return ml_make_nan();
    if (ML_UNLIKELY(!ml_isfinite(s) || !ml_isfinite(a))) return ml_make_nan();
    if (!(s > a)) return ml_make_nan();
    return 1.0 / (s - a);
}

ML_API double ml_laplace_sin(double s, double w) {
    /* L{sin(wt)} = w/(s^2+w^2), s>0. */
    if (ML_UNLIKELY(ml_isnan(s) || ml_isnan(w))) return ml_make_nan();
    if (ML_UNLIKELY(!ml_isfinite(s) || !ml_isfinite(w))) return ml_make_nan();
    if (!(s > 0.0)) return ml_make_nan();
    {
        long double den = (long double)s * (long double)s + (long double)w * (long double)w;
        if (den == 0.0L) return ml_make_nan();
        double r = (double)((long double)w / den);
        return ml_isfinite(r) ? r : ml_make_nan();
    }
}

ML_API double ml_laplace_cos(double s, double w) {
    /* L{cos(wt)} = s/(s^2+w^2), s>0. */
    if (ML_UNLIKELY(ml_isnan(s) || ml_isnan(w))) return ml_make_nan();
    if (ML_UNLIKELY(!ml_isfinite(s) || !ml_isfinite(w))) return ml_make_nan();
    if (!(s > 0.0)) return ml_make_nan();
    {
        long double den = (long double)s * (long double)s + (long double)w * (long double)w;
        if (den == 0.0L) return ml_make_nan();
        double r = (double)((long double)s / den);
        return ml_isfinite(r) ? r : ml_make_nan();
    }
}
