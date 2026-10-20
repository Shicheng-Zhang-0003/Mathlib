#include "ml_interval.h"
#include "ml_core.h"

#include <float.h>
#include <fenv.h>
#include <math.h>
int ml_interval_is_supported (void) {
#if defined(FE_TONEAREST)
    return FLT_RADIX == 2 && DBL_MANT_DIG == 53 && DBL_MAX_EXP == 1024 &&
           fegetround () == FE_TONEAREST;
#else
    return 0;
#endif
}
int ml_interval_is_valid (ml_interval_t x) {
    return ml_isfinite (x.lo) && ml_isfinite (x.hi) && x.lo <= x.hi;
}
ml_status_t ml_interval_make (double lo, double hi, ml_interval_t *out) {
    ml_interval_t r;
    if (out == NULL || !ml_isfinite (lo) || !ml_isfinite (hi) || lo > hi)
        return ML_ERR_INVALID_ARG;
    r.lo = lo;
    r.hi = hi;
    *out = r;
    return ML_SUCCESS;
}
ml_status_t ml_interval_point (double x, ml_interval_t *out) {
    return ml_interval_make (x, x, out);
}
/* Widen a nearest-rounded result by one adjacent binary64 value in each
 * direction. If either widened endpoint is infinite, this finite-endpoint
 * representation cannot hold the enclosure, so fail without writing out. */
static ml_status_t iv_widen (double lo, double hi, ml_interval_t *out) {
    ml_interval_t r;
    double lower, upper;
    if (!ml_isfinite (lo) || !ml_isfinite (hi)) return ML_ERR_OVERFLOW;
    lower = nextafter (lo, -ml_make_inf (0));
    upper = nextafter (hi, ml_make_inf (0));
    if (!ml_isfinite (lower) || !ml_isfinite (upper)) return ML_ERR_OVERFLOW;
    r.lo = lower;
    r.hi = upper;
    *out = r;
    return ML_SUCCESS;
}
static int iv_inputs_ok (ml_interval_t a, ml_interval_t b,
                         ml_interval_t *out) {
    return out != NULL && ml_interval_is_valid (a) && ml_interval_is_valid (b);
}
ml_status_t ml_interval_add (ml_interval_t a, ml_interval_t b,
                             ml_interval_t *out) {
    if (!iv_inputs_ok (a, b, out)) return ML_ERR_INVALID_ARG;
    if (!ml_interval_is_supported ()) return ML_ERR_UNSUPPORTED;
    return iv_widen (a.lo + b.lo, a.hi + b.hi, out);
}
ml_status_t ml_interval_sub (ml_interval_t a, ml_interval_t b,
                             ml_interval_t *out) {
    if (!iv_inputs_ok (a, b, out)) return ML_ERR_INVALID_ARG;
    if (!ml_interval_is_supported ()) return ML_ERR_UNSUPPORTED;
    return iv_widen (a.lo - b.hi, a.hi - b.lo, out);
}
ml_status_t ml_interval_mul (ml_interval_t a, ml_interval_t b,
                             ml_interval_t *out) {
    double p0, p1, p2, p3, lo, hi;
    if (!iv_inputs_ok (a, b, out)) return ML_ERR_INVALID_ARG;
    if (!ml_interval_is_supported ()) return ML_ERR_UNSUPPORTED;
    p0 = a.lo * b.lo;
    p1 = a.lo * b.hi;
    p2 = a.hi * b.lo;
    p3 = a.hi * b.hi;
    if (!ml_isfinite (p0) || !ml_isfinite (p1) ||
        !ml_isfinite (p2) || !ml_isfinite (p3)) return ML_ERR_OVERFLOW;
    lo = p0;
    if (p1 < lo) lo = p1;
    if (p2 < lo) lo = p2;
    if (p3 < lo) lo = p3;
    hi = p0;
    if (p1 > hi) hi = p1;
    if (p2 > hi) hi = p2;
    if (p3 > hi) hi = p3;
    return iv_widen (lo, hi, out);
}
ml_status_t ml_interval_div (ml_interval_t a, ml_interval_t b,
                             ml_interval_t *out) {
    ml_interval_t reciprocal, result;
    ml_status_t st;
    if (!iv_inputs_ok (a, b, out)) return ML_ERR_INVALID_ARG;
    if (!ml_interval_is_supported ()) return ML_ERR_UNSUPPORTED;
    if (b.lo <= 0.0 && b.hi >= 0.0) return ML_ERR_SINGULAR;
    st = iv_widen (1.0 / b.hi, 1.0 / b.lo, &reciprocal);
    if (st != ML_SUCCESS) return st;
    st = ml_interval_mul (a, reciprocal, &result);
    if (st != ML_SUCCESS) return st;
    *out = result;
    return ML_SUCCESS;
}
ml_status_t ml_interval_neg (ml_interval_t a, ml_interval_t *out) {
    ml_interval_t r;
    if (out == NULL || !ml_interval_is_valid (a)) return ML_ERR_INVALID_ARG;
    r.lo = -a.hi;
    r.hi = -a.lo;
    *out = r;
    return ML_SUCCESS;
}
ml_status_t ml_interval_hull (ml_interval_t a, ml_interval_t b,
                              ml_interval_t *out) {
    ml_interval_t r;
    if (!iv_inputs_ok (a, b, out)) return ML_ERR_INVALID_ARG;
    r.lo = a.lo < b.lo ? a.lo : b.lo;
    r.hi = a.hi > b.hi ? a.hi : b.hi;
    *out = r;
    return ML_SUCCESS;
}
int ml_interval_contains (ml_interval_t x, double value) {
    return ml_interval_is_valid (x) && !ml_isnan (value) &&
           x.lo <= value && value <= x.hi;
}
int ml_interval_contains_zero (ml_interval_t x) {
    return ml_interval_is_valid (x) && x.lo <= 0.0 && x.hi >= 0.0;
}
