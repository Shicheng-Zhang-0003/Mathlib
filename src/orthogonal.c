#include "ml_compiler.h"
#include "ml_orthogonal.h"

ML_API double ml_legendre_p(int n, double x) {
    if (ml_isnan(x) || ml_isinf(x)) return ml_make_nan();
    if (n < 0) return ml_make_nan();
    if (n == 0) return 1.0;
    if (n == 1) return x;
    if (x == 1.0) return 1.0;
    if (x == -1.0) return (n % 2 == 0) ? 1.0 : -1.0;
    {
        long double xl = (long double)x;
        long double p0 = 1.0L, p1 = xl, p2 = 0.0L;
        for (int k = 1; k < n; k++) {
            /* (k+1)P_{k+1} = (2k+1)xP_k - kP_{k-1} */
            p2 = ((2.0L * (long double)k + 1.0L) * xl * p1
                  - (long double)k * p0) / (long double)(k + 1);
            p0 = p1; p1 = p2;
        }
        double r = (double)p1;
        return ml_isfinite(r) ? r : ml_make_nan();
    }
}

ML_API double ml_legendre_p_deriv(int n, double x) {
    if (ml_isnan(x) || ml_isinf(x)) return ml_make_nan();
    if (n < 0) return ml_make_nan();
    if (n == 0) return 0.0;
    if (x == 1.0) return 0.5 * (double)n * (double)(n + 1);
    if (x == -1.0) {
        double v = 0.5 * (double)n * (double)(n + 1);
        return (n % 2 == 0) ? -v : v;
    }
    {
        /* (x^2-1)P'_n = n(xP_n - P_{n-1}) */
        long double xl = (long double)x;
        long double p0 = 1.0L, p1 = xl, p2 = 0.0L;
        if (n == 1) return 1.0;
        for (int k = 1; k < n; k++) {
            p2 = ((2.0L * (long double)k + 1.0L) * xl * p1
                  - (long double)k * p0) / (long double)(k + 1);
            p0 = p1; p1 = p2;
        }
        long double d = (long double)n * (xl * p1 - p0) / (xl * xl - 1.0L);
        double r = (double)d;
        return ml_isfinite(r) ? r : ml_make_nan();
    }
}

ML_API double ml_chebyshev_t(int n, double x) {
    if (ml_isnan(x) || ml_isinf(x)) return ml_make_nan();
    if (n < 0) return ml_make_nan();
    if (n == 0) return 1.0;
    if (n == 1) return x;
    {
        long double xl = (long double)x;
        long double t0 = 1.0L, t1 = xl, t2 = 0.0L;
        for (int k = 1; k < n; k++) {
            t2 = 2.0L * xl * t1 - t0;
            t0 = t1; t1 = t2;
        }
        double r = (double)t1;
        return ml_isfinite(r) ? r : ml_make_nan();
    }
}

ML_API double ml_chebyshev_u(int n, double x) {
    if (ml_isnan(x) || ml_isinf(x)) return ml_make_nan();
    if (n < 0) return ml_make_nan();
    if (n == 0) return 1.0;
    if (n == 1) return 2.0 * x;
    {
        long double xl = (long double)x;
        long double u0 = 1.0L, u1 = 2.0L * xl, u2 = 0.0L;
        for (int k = 1; k < n; k++) {
            u2 = 2.0L * xl * u1 - u0;
            u0 = u1; u1 = u2;
        }
        double r = (double)u1;
        return ml_isfinite(r) ? r : ml_make_nan();
    }
}

ML_API double ml_hermite_h(int n, double x) {
    if (ml_isnan(x) || ml_isinf(x)) return ml_make_nan();
    if (n < 0) return ml_make_nan();
    if (n == 0) return 1.0;
    if (n == 1) return 2.0 * x;
    {
        long double xl = (long double)x;
        long double h0 = 1.0L, h1 = 2.0L * xl, h2 = 0.0L;
        for (int k = 1; k < n; k++) {
            h2 = 2.0L * xl * h1 - 2.0L * (long double)k * h0;
            h0 = h1; h1 = h2;
        }
        double r = (double)h1;
        return ml_isfinite(r) ? r : ml_make_nan();
    }
}

ML_API double ml_laguerre_l(int n, double x) {
    if (ml_isnan(x) || ml_isinf(x)) return ml_make_nan();
    if (n < 0) return ml_make_nan();
    if (x < 0.0) return ml_make_nan();
    if (n == 0) return 1.0;
    if (n == 1) return 1.0 - x;
    {
        long double xl = (long double)x;
        long double l0 = 1.0L, l1 = 1.0L - xl, l2 = 0.0L;
        for (int k = 1; k < n; k++) {
            l2 = ((2.0L * (long double)k + 1.0L - xl) * l1
                  - (long double)k * l0) / (long double)(k + 1);
            l0 = l1; l1 = l2;
        }
        double r = (double)l1;
        return ml_isfinite(r) ? r : ml_make_nan();
    }
}

ML_API double ml_laguerre_l_gen(int n, double alpha, double x) {
    if (ml_isnan(x) || ml_isnan(alpha) || ml_isinf(x) || ml_isinf(alpha)) return ml_make_nan();
    if (n < 0 || !(alpha > -1.0) || !ml_isfinite(alpha)) return ml_make_nan();
    if (x < 0.0) return ml_make_nan();
    if (n == 0) return 1.0;
    if (n == 1) return 1.0 + alpha - x;
    {
        long double xl = (long double)x, al = (long double)alpha;
        long double l0 = 1.0L, l1 = 1.0L + al - xl, l2 = 0.0L;
        for (int k = 1; k < n; k++) {
            l2 = ((2.0L * (long double)k + 1.0L + al - xl) * l1
                  - ((long double)k + al) * l0) / (long double)(k + 1);
            l0 = l1; l1 = l2;
        }
        double r = (double)l1;
        return ml_isfinite(r) ? r : ml_make_nan();
    }
}
