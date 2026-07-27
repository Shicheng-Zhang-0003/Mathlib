#ifndef LIBMATHC_IEEE754_H
#define LIBMATHC_IEEE754_H

#include <stdint.h>
#include <string.h>
#include "ml_core.h"

/* ============================================================================
 * v11S AUDIT IP-1:
 * These helpers are still experimental, but they are now safe against:
 *   - unbounded integer casts
 *   - NaN / Inf misuse
 *   - obvious mantissa-balancing bugs
 * ========================================================================== */

ML_INLINE double logarithm_ieee754(double x) {
    if (ml_isnan(x)) return x;
    if (x == 0.0) return -ml_make_inf(0);
    if (x < 0.0) return ml_make_nan();
    if (ml_isinf(x)) return x;
    if (x == 1.0) return 0.0;

    uint64_t bits;
    memcpy(&bits, &x, sizeof(uint64_t));

    uint64_t exp_bits = (bits >> 52) & 0x7FFULL;
    uint64_t mantissa = bits & 0x000FFFFFFFFFFFFFULL;

    int e;
    double m;

    if (exp_bits == 0) {
        /* Subnormal */
        if (mantissa == 0) return -ml_make_inf(0);

        m = (double)mantissa / 4503599627370496.0; /* 2^52 */
        e = -1022;

        while (m < 1.0) {
            m *= 2.0;
            e--;
        }
    } else {
        mantissa |= 0x10000000000000ULL;
        m = (double)mantissa / 4503599627370496.0; /* 2^52 */
        e = (int)exp_bits - 1023;
    }

    /* Balance into [sqrt(1/2), sqrt(2)] */
    if (m > 1.4142135623730951) {
        m *= 0.5;
        e++;
    }

    double z = (m - 1.0) / (m + 1.0);
    double z2 = z * z;

    double res = z;
    double term = z;

    for (int i = 3; i <= 15; i += 2) {
        term *= z2;
        res += term / i;
    }

    return 2.0 * res + (double)e * ML_LN2;
}

ML_INLINE double exponential_ieee754(double x) {
    if (ml_isnan(x)) return x;
    if (x == 0.0) return 1.0;
    if (ml_isinf(x)) return x > 0.0 ? ml_make_inf(0) : 0.0;

    if (x > 709.782712893384) return ml_make_inf(0);
    if (x < -745.133219101941) return 0.0;

    double n_d = x / ML_LN2;
    int n = (int)(n_d + (n_d > 0.0 ? 0.5 : -0.5));

    double r = x - (double)n * ML_LN2;

    double res = 1.0;
    double term = 1.0;

    for (int i = 1; i <= 15; i++) {
        term *= r / i;
        res += term;
    }

    uint64_t bits;
    memcpy(&bits, &res, sizeof(uint64_t));

    int64_t current_exp = (int64_t)((bits >> 52) & 0x7FFULL);
    int64_t new_exp = current_exp + n;

    if (new_exp <= 0) return 0.0;
    if (new_exp >= 0x7FF) return ml_make_inf(0);

    bits = (bits & 0x800FFFFFFFFFFFFFULL) | ((uint64_t)new_exp << 52);

    double final_res;
    memcpy(&final_res, &bits, sizeof(double));

    return final_res;
}

#endif /* LIBMATHC_IEEE754_H */
