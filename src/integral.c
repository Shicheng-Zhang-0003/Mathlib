#include "ml_compiler.h"
#include "ml_integral.h"
#include "ml_trig.h"

/* MATHLIB_V12A1_GAMMA_FIX_RESTORE */
/*
 * v11S contract:
 * ml_factorial_float(x) means x! = Gamma(x + 1).
 */
ML_API double ml_factorial_float(double x) {
    if (ml_isnan(x)) return x;
    if (x < 0.0) return ml_make_nan();
    if (ml_isinf(x)) return ml_make_inf(0);
    if (x == 0.0) return 1.0;
    return ml_gamma_new(x + 1.0);
}

/*
 * Educational / experimental numerical integrator.
 * Retained for compatibility.
 */
ML_API double ml_integral_traditional(double a, double b, double exponent, double additive, double d) {
    if (ml_isnan(a) || ml_isnan(b) || ml_isnan(exponent) ||
        ml_isnan(additive) || ml_isnan(d)) {
        return ml_make_nan();
    }
    if (ml_isinf(a) || ml_isinf(b) || ml_isinf(exponent) ||
        ml_isinf(additive) || ml_isinf(d)) {
        return ml_make_nan();
    }
    if (d == 0.0) {
        return ml_make_nan();
    }
    if ((d > 0.0 && a >= b) || (d < 0.0 && a <= b)) {
        return 0.0;
    }
    double result = 0.0;
    double x = a;
    const int max_steps = 10000000;
    for (int step = 0; step < max_steps; step++) {
        if ((d > 0.0 && x >= b) || (d < 0.0 && x <= b)) {
            return result;
        }
        double term = ml_pow(x, exponent) + additive;
        if (ml_isnan(term)) {
            return ml_make_nan();
        }
        result += term * d;
        double next_x = x + d;
        if (next_x == x) {
            return ml_make_nan();
        }
        x = next_x;
    }
    return ml_make_nan();
}

/* MATHLIB_V12A1_GAMMA_LANCZOS */
/*
 * Lanczos approximation for the gamma function.
 *
 * Uses g=7, n=9 coefficients (Numerical Recipes / Godfrey).
 * Achieves ~15 significant digits across the positive real line.
 *
 * The old v11S implementation used a degree-8 polynomial on [1,2]
 * with iterative reduction, giving only ~1e-2 relative accuracy.
 * This is a 13-order-of-magnitude improvement.
 */

static const double ml_lanczos_coeff[9] = {
     0.99999999999980993,
   676.5203681218851,
  -1259.1392167224028,
   771.32342877765313,
  -176.61502916214059,
    12.507343278686905,
    -0.13857109526572012,
     9.9843695780195716e-6,
     1.5056327351493116e-7
};

#define ML_LANCZOS_G 7.0
#define ML_GAMMA_OVERFLOW 171.6243769563027

/*
 * Internal: log Gamma(x) via Lanczos for x > 0.
 *
 * log Gamma(x) = 0.5*log(2pi) + (z+0.5)*log(t) - t + log(Ag(z))
 * where z = x - 1, t = z + g + 0.5.
 */
static double ml_lgamma_positive(double x) {
    /* MATHLIB_V12A1_GAMMA_KAHAN */
    double z = x - 1.0;

    /* Kahan summation for the Lanczos coefficient sum */
    double ag = ml_lanczos_coeff[0];
    double comp = 0.0;
    for (int i = 1; i < 9; i++) {
        double term = ml_lanczos_coeff[i] / (z + (double)i);
        double y_k = term - comp;
        double t_k = ag + y_k;
        comp = (t_k - ag) - y_k;
        ag = t_k;
    }

    double t = z + ML_LANCZOS_G + 0.5;
        /* MATHLIB_V12A1_GAMMA_LOG_SPLIT */
    /*
     * Double-double log to avoid (z+0.5)*log(t) amplification.
     *
     * The old code computed (z+0.5) * ml_log(t), where ml_log(t)
     * has ~1 ULP error. When multiplied by (z+0.5) which can be
     * up to ~171, the error is amplified to ~(z+0.5) ULP.
     *
     * Fix: compute log(t) as double-double (log_hi + log_lo),
     * then multiply by (z+0.5) using FMA for extended precision.
     */
    double log_hi, log_lo;
    ml_log_split(t, &log_hi, &log_lo);
    double zp5 = z + 0.5;
    double prod_hi = zp5 * log_hi;
    double prod_lo = ML_FMA(zp5, log_hi, -prod_hi) + zp5 * log_lo;
    double prod = prod_hi + prod_lo;
    return 0.91893853320467274178 + prod - t + ml_log(ag);
}

ML_API double ml_lgamma(double x) {
    /* MATHLIB_V12A1_LGAMMA */
    if (ml_isnan(x)) return x;
    if (ml_isinf(x)) return ml_make_inf(0);

    /* Poles at zero and negative integers */
    if (x <= 0.0 && x == ml_round(x)) return ml_make_inf(0);

    if (x > 0.0) {
        return ml_lgamma_positive(x);
    }

    /* Reflection: log|Gamma(x)| = log(pi) - log|sin(pi*x)| - log Gamma(1-x) */
    double sinpx = ml_sin(ML_PI * x);
    if (sinpx == 0.0) return ml_make_inf(0);
    return ml_log(ML_PI / ml_fabs(sinpx)) - ml_lgamma_positive(1.0 - x);
}

ML_API double ml_gamma_new(double x) {
    /* MATHLIB_V12A1_GAMMA_LANCZOS */
    if (ml_isnan(x)) return x;
    if (ml_isinf(x)) return x > 0.0 ? ml_make_inf(0) : ml_make_nan();

    /* Poles at zero and negative integers */
    if (x <= 0.0 && x == ml_round(x)) return ml_make_nan();

    if (x > 0.0) {
        if (x > ML_GAMMA_OVERFLOW) return ml_make_inf(0);
        /* MATHLIB_V12A1_GAMMA_DIRECT */
        /* Direct computation: avoids log(ag) -> exp roundtrip */
        {
            double z = x - 1.0;
            double ag = ml_lanczos_coeff[0];
            double comp = 0.0;
            for (int i = 1; i < 9; i++) {
                double term = ml_lanczos_coeff[i] / (z + (double)i);
                double y_k = term - comp;
                double t_k = ag + y_k;
                comp = (t_k - ag) - y_k;
                ag = t_k;
            }
            double t = z + ML_LANCZOS_G + 0.5;
            /* MATHLIB_V12A1_GAMMA_LOG_SPLIT */
            /*
             * Double-double log to avoid (z+0.5)*log(t) amplification.
             */
            double log_hi, log_lo;
            ml_log_split(t, &log_hi, &log_lo);
            double zp5 = z + 0.5;
            double prod_hi = zp5 * log_hi;
            double prod_lo = ML_FMA(zp5, log_hi, -prod_hi) + zp5 * log_lo;
            double log_part = (prod_hi + prod_lo) - t;
            /* sqrt(2*pi) = 2.50662827463100050242 */
            return 2.50662827463100050242 * ag * ml_exp(log_part);
        }
    }

    /* Reflection: Gamma(x) = pi / (sin(pi*x) * Gamma(1-x)) */
    double sinpx = ml_sin(ML_PI * x);
    if (sinpx == 0.0) return ml_make_nan();
    if (1.0 - x > ML_GAMMA_OVERFLOW) {
        /* Gamma(1-x) overflows, so Gamma(x) -> 0 */
        return ml_copysign(0.0, sinpx);
    }
    return ML_PI / (sinpx * ml_exp(ml_lgamma_positive(1.0 - x)));
}
