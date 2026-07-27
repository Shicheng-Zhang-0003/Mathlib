#include "ml_compiler.h"
#include "ml_trig.h"
#include "internal/minimax.h"
#include "internal/cordic.h"
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
        return (s > 0) ? (ml_make_inf(0)) : (-ml_make_inf(0));
    }
    return s / c;
}

ML_API double ml_atan(double x) {
    if (ml_isnan(x)) return x;
    if (ml_isinf(x)) return ml_copysign(ML_PI / 2.0, x);

    if (x > 1.0) return (ML_PI / 2.0) - ml_atan(1.0 / x);
    if (x < -1.0) return -(ML_PI / 2.0) - ml_atan(1.0 / x);
    if (x > 0.5) return (ML_PI / 4.0) + ml_atan((x - 1.0) / (x + 1.0));
    if (x < -0.5) return -(ML_PI / 4.0) + ml_atan((x + 1.0) / (1.0 - x));

    double result = x, term = x, x2 = x * x;
    for (int i = 3; i <= 21; i += 2) { term *= -x2; result += term / i; }
    return result;
}

ML_API double ml_asin(double x) {
    if (x < -1.0 || x > 1.0) return ml_make_nan();
    return 2.0 * ml_atan(x / (1.0 + ml_sqrt(1.0 - x * x)));
}

ML_API double ml_acos(double x) {
    if (x < -1.0 || x > 1.0) return ml_make_nan();
    return (ML_PI / 2.0) - ml_asin(x);
}

ML_API double ml_acot(double x) {
    return (ML_PI / 2.0) - ml_atan(x);
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
