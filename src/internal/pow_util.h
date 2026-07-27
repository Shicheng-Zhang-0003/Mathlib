#ifndef ML_INTERNAL_POW_UTIL_H
#define ML_INTERNAL_POW_UTIL_H

#include "ml_core.h"

/*
 * MATHLIB_CLOSURE_P0:
 * Shared exponent-classification helpers for real and complex power paths.
 */

static inline int ml_is_integer_double(double y) {
    if (ml_isnan(y) || ml_isinf(y)) return 0;
    if (y == 0.0) return 1;

    /*
     * At or beyond 2^53, all representable finite doubles are integral.
     */
    if (ml_fabs(y) >= 9007199254740992.0) return 1;

    return ml_round(y) == y;
}

static inline int ml_is_odd_integer_double(double y) {
    long long yi;
    long long ay;

    if (!ml_is_integer_double(y)) return 0;

    /*
     * At or beyond 2^53, representable integers are multiples of two
     * or larger powers of two, so they cannot be odd.
     */
    if (ml_fabs(y) >= 9007199254740992.0) return 0;

    yi = (long long)y;
    ay = (yi < 0) ? -yi : yi;

    return (int)(ay % 2LL);
}

#endif /* ML_INTERNAL_POW_UTIL_H */
