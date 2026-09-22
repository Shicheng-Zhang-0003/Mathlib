#include "ml_compiler.h"
#include "ml_trig.h"
#include "internal/minimax.h"
#include "internal/cordic.h"
#include "internal/hypot.h"
#include "ml_fixed_point.h"

ML_API double ml_sin(double x) {
    if (ml_isnan(x) || ml_isinf(x)) return ml_make_nan();
#if defined(MATHLIB_PROFILE_EMBEDDED)
    x = ml_fmod(x, 2.0 * ML_PI);
    ml_q16_16_t f_in = (ml_q16_16_t)(x * 65536.0);
    ml_q16_16_t s, c;
    ml_cordic_sincos_fixed(f_in, &s, &c);
    return (double)s / 65536.0;
#else
    return ml_minimax_sin(x);
#endif
}

ML_API double ml_cos(double x) {
    if (ml_isnan(x) || ml_isinf(x)) return ml_make_nan();
#if defined(MATHLIB_PROFILE_EMBEDDED)
    x = ml_fmod(x, 2.0 * ML_PI);
    ml_q16_16_t f_in = (ml_q16_16_t)(x * 65536.0);
    ml_q16_16_t s, c;
    ml_cordic_sincos_fixed(f_in, &s, &c);
    return (double)c / 65536.0;
#else
    return ml_minimax_cos(x);
#endif
}

ML_API double ml_tan(double x) {
    double s = ml_sin(x);
    double c = ml_cos(x);
    if (c == 0.0) {
        if (s == 0.0) return ml_make_nan(); /* NaN for 0/0 */
        /* IEEE s/c gives correctly-signed Inf from signed zeros. */
        return s / c;
    }
    return s / c;
}

#ifndef ML_PI_HI_D
#define ML_PI_HI_D 0x1.921fb54442d18p+1
#endif
#ifndef ML_PI_LO_D
#define ML_PI_LO_D 0x1.1a62633145c07p-53
#endif

ML_API double ml_sinpi(double x) {
    if (ml_isnan(x) || ml_isinf(x)) return ml_make_nan();
    if (x == 0.0) return x;
    /* Reduce mod 2: sin(pi*x) period 2. Use fmod for exact reduction,
     * then compensated pi*r = r*HI + r*LO. Handles large |x| via fmod. */
    {
        double r = ml_fmod(x, 2.0);
        if (r > 1.0) r -= 2.0;
        if (r < -1.0) r += 2.0;
        if (r == 0.0) return ml_copysign(0.0, x);
        if (r == 1.0 || r == -1.0) return ml_copysign(0.0, x);
        if (r == 0.5) return 1.0;
        if (r == -0.5) return -1.0;
        {
            double ahi = r * ML_PI_HI_D;
            double alo = ML_FMA(r, ML_PI_HI_D, -ahi) + r * ML_PI_LO_D;
            double arg = ahi + alo;
            if (!ml_isfinite(arg)) return ml_sin(ML_PI * x);
            return ml_sin(arg);
        }
    }
}

ML_API double ml_cospi(double x) {
    if (ml_isnan(x) || ml_isinf(x)) return ml_make_nan();
    /* cos(pi*x) = sin(pi*(x+0.5)) with sign care for large values. */
    {
        double ax = ml_fabs(x);
        if (ax == 0.0) return 1.0;
        /* Reduce mod 2 first for stability. */
        double r = ml_fmod(x, 2.0);
        if (r > 1.0) r -= 2.0;
        if (r < -1.0) r += 2.0;
        if (r == 0.5 || r == -0.5) return 0.0;
        {
            double ahi = r * ML_PI_HI_D;
            double alo = ML_FMA(r, ML_PI_HI_D, -ahi) + r * ML_PI_LO_D;
            double arg = ahi + alo;
            if (!ml_isfinite(arg)) return ml_cos(ML_PI * x);
            return ml_cos(arg);
        }
    }
}

static inline ml_ddx_t ml_atan_raw_dd(double x) {
    /* atan(x)=x*P(x^2), P=Σ(-1)^k t^k/(2k+1), 22 terms, DD Horner. */
    double t = x * x;
    ml_ddx_t acc = ml_ddx_from_d(1.0 / 43.0);
    for (int k = 20; k >= 0; k--) {
        acc = ml_ddx_mul_d(acc, t);
        acc = ml_ddx_add_d(acc, ((k % 2 == 0) ? 1.0 : -1.0) / (double)(2 * k + 1));
    }
    return ml_ddx_mul_d(acc, x);
}

ML_API double ml_atan(double x) {
    if (ml_isnan(x)) return x;
    if (ml_isinf(x)) return ml_copysign(ML_PI / 2.0, x);

    /* DD pi slices (HI/LO exact halves/quarters). */
    {
        double pih = ML_PI_HI_D, pil = ML_PI_LO_D;
        double pi2h = pih * 0.5, pi2l = pil * 0.5;
        double pi4h = pih * 0.25, pi4l = pil * 0.25;
        if (x > 1.0) {
            double a = ml_atan(1.0 / x);
            ml_ddx_t r = ml_ddx_from_d(a);
            r = ml_ddx_renorm(pi2h - r.hi, pi2l - r.lo);
            /* pi/2 - a via DD: hi=pi2h-a with correction. Use generic: */
            {
                double s, e;
                s = ml_two_sum(pi2h, -a, &e);
                e += pi2l;
                return ml_ddx_to_d(ml_ddx_renorm(s, e));
            }
        }
        if (x < -1.0) {
            double a = ml_atan(1.0 / x);
            double s, e;
            s = ml_two_sum(-pi2h, -a, &e);
            e += -pi2l;
            return ml_ddx_to_d(ml_ddx_renorm(s, e));
        }
        if (x > 0.5) {
            double a = ml_atan((x - 1.0) / (x + 1.0));
            double s, e;
            s = ml_two_sum(pi4h, a, &e);
            e += pi4l;
            return ml_ddx_to_d(ml_ddx_renorm(s, e));
        }
        if (x < -0.5) {
            double a = ml_atan((x + 1.0) / (1.0 - x));
            double s, e;
            s = ml_two_sum(-pi4h, a, &e);
            e += -pi4l;
            return ml_ddx_to_d(ml_ddx_renorm(s, e));
        }
        if (x > 0.2679491924311227) {
            double t = 0.2679491924311227;
            double a = ml_atan((x - t) / (1.0 + t * x));
            /* pi/12 in DD to match pi/2, pi/4 paths (was plain double). */
            double c12h = ML_PI_HI_D / 12.0;
            double c12l = ML_PI_LO_D / 12.0;
            double s, e;
            s = ml_two_sum(c12h, a, &e);
            e += c12l;
            return ml_ddx_to_d(ml_ddx_renorm(s, e));
        }
        if (x < -0.2679491924311227) {
            double t = 0.2679491924311227;
            double a = ml_atan((x + t) / (1.0 - t * x));
            double c12h = ML_PI_HI_D / 12.0;
            double c12l = ML_PI_LO_D / 12.0;
            double s, e;
            s = ml_two_sum(-c12h, a, &e);
            e += -c12l;
            return ml_ddx_to_d(ml_ddx_renorm(s, e));
        }
    }

    return ml_ddx_to_d(ml_atan_raw_dd(x));
}

ML_API double ml_asin(double x) {
    if (x < -1.0 || x > 1.0) return ml_make_nan();
    if (x == 1.0) return ML_PI / 2.0;
    if (x == -1.0) return -ML_PI / 2.0;
    if (x > 0.9) return (ML_PI / 2.0) - 2.0 * ml_atan(ml_sqrt((1.0 - x) / (1.0 + x)));
    if (x < -0.9) return -(ML_PI / 2.0) + 2.0 * ml_atan(ml_sqrt((1.0 + x) / (1.0 - x)));
    return 2.0 * ml_atan(x / (1.0 + ml_sqrt(1.0 - x * x)));
}

ML_API double ml_acos(double x) {
    if (x < -1.0 || x > 1.0) return ml_make_nan();
    if (x >= 0.0) return 2.0 * ml_atan(ml_sqrt((1.0 - x) / (1.0 + x)));
    return ML_PI - 2.0 * ml_atan(ml_sqrt((1.0 + x) / (1.0 - x)));
}

ML_API double ml_acot(double x) {
    if (ml_isnan(x)) return x;
    if (ml_isinf(x)) return (x > 0.0) ? 0.0 : ML_PI;
    if (x == 0.0) return ML_PI / 2.0;
    /* For |x|>1 use atan(1/x) form to avoid catastrophic cancellation
     * in pi/2-atan(x) (atan(x)~pi/2-1/x loses relative accuracy).
     * acot(x)=atan(1/x) for x>0, pi+atan(1/x) for x<0. */
    if (x > 1.0) return ml_atan(1.0 / x);
    if (x < -1.0) return ML_PI + ml_atan(1.0 / x);
    return (ML_PI / 2.0) - ml_atan(x);
}

ML_API double ml_sec(double x) {
    double c = ml_cos(x);
    double s = ml_sin(x);
    if (c == 0.0) {
        /* 0/0 (both zero) only from rounding failure; else IEEE 1/c
         * gives correctly-signed Inf from signed zero. */
        if (s == 0.0) return ml_make_nan();
        return 1.0 / c;
    }
    return 1.0 / c;
}

ML_API double ml_csc(double x) {
    double s = ml_sin(x);
    double c = ml_cos(x);
    if (s == 0.0) {
        if (c == 0.0) return ml_make_nan();
        return 1.0 / s;
    }
    return 1.0 / s;
}

ML_API double ml_cot(double x) {
    double s = ml_sin(x);
    double c = ml_cos(x);
    if (s == 0.0) {
        if (c == 0.0) return ml_make_nan();
        /* IEEE c/s preserves signs of both zeros. */
    }
    return c / s;
}

ML_API double ml_sinc(double x) {
    if (ml_isnan(x)) return x;
    if (ml_isinf(x)) return 0.0;
    if (x == 0.0) return 1.0;
    {
        double ax = ml_fabs(x);
        if (ax < 1.5e-8) return 1.0 - x * x / 6.0;
        return ml_sin(x) / x;
    }
}

/* FIX: Strict IEEE-754 Annex F compliance for signed zeros and Infinities.
 * ml_atan2(+0.0, -0.0) now correctly returns PI instead of 0.0. */
ML_API double ml_atan2(double y, double x) {
    /* v11S AUDIT IP-2: NaN propagation + signed-zero / infinity fixes */
    if (ml_isnan(x) || ml_isnan(y)) return ml_make_nan();

    int y_neg = ml_signbit(y);
    int x_neg = ml_signbit(x);
    int y_zero = (y == 0.0);
    int x_zero = (x == 0.0);

    if (y_zero && x_zero) {
        if (y_neg && x_neg) return -ML_PI;
        if (y_neg && !x_neg) return -0.0;
        if (!y_neg && x_neg) return ML_PI;
        return 0.0;
    }

    if (x_zero) {
        return y_neg ? -ML_PI / 2.0 : ML_PI / 2.0;
    }

    if (y_zero) {
        if (x_neg) return y_neg ? -ML_PI : ML_PI;
        return y_neg ? -0.0 : 0.0;
    }

    if (ml_isinf(x) && ml_isinf(y)) {
        double pi_4 = ML_PI / 4.0;

        if (y_neg && x_neg) return -3.0 * pi_4;
        if (y_neg && !x_neg) return -pi_4;
        if (!y_neg && x_neg) return 3.0 * pi_4;
        return pi_4;
    }

    if (ml_isinf(x)) {
        if (x_neg) return ml_copysign(ML_PI, y);
        return ml_copysign(0.0, y);
    }

    if (ml_isinf(y)) {
        return y_neg ? -ML_PI / 2.0 : ML_PI / 2.0;
    }

    double a = ml_atan(y / x);

    if (x_neg) {
        return y_neg ? a - ML_PI : a + ML_PI;
    }

    return a;
}

ML_API double ml_r_form(double a, double b, double *R, double *alpha) {
    if (ml_isnan(a) || ml_isnan(b)) return ml_make_nan();
    if (!ml_isfinite(a) || !ml_isfinite(b)) return ml_make_nan();
    {
        double r = ml_hypot_internal(a, b);
        double al = ml_atan2(b, a);
        if (R) *R = r;
        if (alpha) *alpha = al;
        return r;
    }
}

ML_API double ml_triangle_area_sas(double a, double b, double C) {
    if (ml_isnan(a) || ml_isnan(b) || ml_isnan(C)) return ml_make_nan();
    if (!ml_isfinite(a) || !ml_isfinite(b) || !ml_isfinite(C)) return ml_make_nan();
    if (a < 0.0 || b < 0.0) return ml_make_nan();
    return 0.5 * a * b * ml_sin(C);
}

ML_API double ml_law_cos_side(double a, double b, double C) {
    if (ml_isnan(a) || ml_isnan(b) || ml_isnan(C)) return ml_make_nan();
    if (!ml_isfinite(a) || !ml_isfinite(b) || !ml_isfinite(C)) return ml_make_nan();
    if (a < 0.0 || b < 0.0) return ml_make_nan();
    {
        double c2 = a * a + b * b - 2.0 * a * b * ml_cos(C);
        if (c2 < 0.0) c2 = 0.0;
        return ml_sqrt(c2);
    }
}

ML_API double ml_law_sin_side(double a, double A, double B) {
    if (ml_isnan(a) || ml_isnan(A) || ml_isnan(B)) return ml_make_nan();
    if (!ml_isfinite(a) || !ml_isfinite(A) || !ml_isfinite(B)) return ml_make_nan();
    if (a < 0.0) return ml_make_nan();
    {
        double sA = ml_sin(A), sB = ml_sin(B);
        if (sA == 0.0 || !ml_isfinite(sA) || !ml_isfinite(sB)) return ml_make_nan();
        return a * sB / sA;
    }
}

ML_API void ml_polar_to_cart(double r, double theta, double *x, double *y) {
    if (!x || !y) return;
    if (ml_isnan(r) || ml_isnan(theta) || !ml_isfinite(r) || !ml_isfinite(theta)) {
        *x = ml_make_nan(); *y = ml_make_nan(); return;
    }
    if (r < 0.0) { *x = ml_make_nan(); *y = ml_make_nan(); return; }
    *x = r * ml_cos(theta);
    *y = r * ml_sin(theta);
}

ML_API void ml_cart_to_polar(double x, double y, double *r, double *theta) {
    if (!r || !theta) return;
    if (ml_isnan(x) || ml_isnan(y) || !ml_isfinite(x) || !ml_isfinite(y)) {
        *r = ml_make_nan(); *theta = ml_make_nan(); return;
    }
    *r = ml_hypot_internal(x, y);
    *theta = ml_atan2(y, x);
}

ML_API double ml_heron(double a, double b, double c) {
    if (ml_isnan(a) || ml_isnan(b) || ml_isnan(c)) return ml_make_nan();
    if (!ml_isfinite(a) || !ml_isfinite(b) || !ml_isfinite(c)) return ml_make_nan();
    if (a <= 0.0 || b <= 0.0 || c <= 0.0) return ml_make_nan();
    if (a + b <= c || b + c <= a || c + a <= b) return ml_make_nan();
    {
        /* Kahan stable: sort x>=y>=z. */
        double x = a, y = b, z = c, t;
        if (x < y) { t = x; x = y; y = t; }
        if (y < z) { t = y; y = z; z = t; }
        if (x < y) { t = x; x = y; y = t; }
        {
            double v = (x + (y + z)) * (z - (x - y)) * (z + (x - y)) * (x + (y - z));
            if (v <= 0.0 || !ml_isfinite(v)) return ml_make_nan();
            return ml_sqrt(v) * 0.25;
        }
    }
}

ML_API double ml_stewart(double a, double b, double c, double m, double n) {
    /* a=m+n opposite A; b,c adjacent; cevian d: b^2 m + c^2 n = a(d^2+mn). */
    if (ml_isnan(a) || ml_isnan(b) || ml_isnan(c) || ml_isnan(m) || ml_isnan(n)) {
        return ml_make_nan();
    }
    if (!ml_isfinite(a) || !ml_isfinite(b) || !ml_isfinite(c) || !ml_isfinite(m) || !ml_isfinite(n)) {
        return ml_make_nan();
    }
    if (a <= 0.0 || b <= 0.0 || c <= 0.0 || m <= 0.0 || n <= 0.0) return ml_make_nan();
    /* Scale-aware cevian closure: relative to max(|m+n|,|a|) so subnormal
     * and huge triangles use the same relative tolerance. */
    {
        double ref = ml_fabs(m + n) > ml_fabs(a) ? ml_fabs(m + n) : ml_fabs(a);
        if (!(ref > 0.0) || !ml_isfinite(ref)) return ml_make_nan();
        if (ml_fabs((m + n) - a) > 1e-12 * ref) return ml_make_nan();
    }
    {
        double d2 = (b * b * m + c * c * n) / a - m * n;
        if (d2 < 0.0) {
            /* Scale-aware clamp: relative to |b^2 m + c^2 n|/a + mn. */
            double scale = ml_fabs((b * b * m + c * c * n) / a) + ml_fabs(m * n);
            if (!(scale > 0.0) || !ml_isfinite(scale)) scale = 1.0;
            if (d2 > -1e-12 * scale) d2 = 0.0;
            else return ml_make_nan();
        }
        return ml_sqrt(d2);
    }
}

ML_API int ml_ceva(double a1, double a2, double b1, double b2, double c1, double c2, double tol) {
    if (ml_isnan(a1) || ml_isnan(a2) || ml_isnan(b1) || ml_isnan(b2) || ml_isnan(c1) || ml_isnan(c2)) {
        return 0;
    }
    if (!ml_isfinite(a1) || !ml_isfinite(a2) || !ml_isfinite(b1) || !ml_isfinite(b2) || !ml_isfinite(c1) || !ml_isfinite(c2)) {
        return 0;
    }
    if (a2 == 0.0 || b2 == 0.0 || c2 == 0.0) return 0;
    if (!(tol > 0.0) || !ml_isfinite(tol)) tol = 1e-9;
    {
        /* Long-double product: double ratios can span 1e-616..1e616,
         * whose product overflows/underflows double; long double
         * (1e4932) holds it exactly for the comparison to 1. */
        long double p = (long double)(a1 / a2) * (long double)(b1 / b2) * (long double)(c1 / c2);
        if (!(p == p)) return 0;
        return (ml_fabs((double)(p - 1.0L)) <= tol) ? 1 : 0;
    }
}
