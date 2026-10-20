#include "ml_compiler.h"
#include "ml_core.h"
#include "ml_math_extra.h"
#include "ml_trig.h"
#include "ml_integral.h"
#include "internal/ieee_exact.h"
#include <limits.h>
#include <math.h>
/* Smallest positive subnormal, constructed exactly. */
static double ml_smallest_subnormal (void) {
    return ml_fp_compose (1ULL, -1074, 0);
}
ML_API double ml_nextafter (double x, double y) {
    uint64_t bits;
    if (ml_isnan (x) || ml_isnan (y)) return ml_make_nan ();
    if (x == y) return y;
    if (x == 0.0) return ml_copysign (ml_smallest_subnormal (), y);

    memcpy (&bits, &x, sizeof (uint64_t));

    /* Step one representable value in the direction of y.  Incrementing or
     * decrementing the raw magnitude bits is exactly the neighbour, because
     * IEEE-754 doubles are monotone in their bit patterns. */
    if (x < y) {
        if (x > 0.0) bits += 1ULL;
        else bits -= 1ULL;
    } else {
        if (x > 0.0) bits -= 1ULL;
        else bits += 1ULL;
    }

    memcpy (&x, &bits, sizeof (double));
    return x;
}
ML_API double ml_nexttoward (double x, long double y) {
    if (ml_isnan (x) || (y != y)) return ml_make_nan ();
    if ((long double) x == y) return x;
    if ((long double) x < y) return ml_nextafter (x, ml_make_inf (0));
    return ml_nextafter (x, ml_make_inf (1));
}
ML_API double ml_modf (double x, double *iptr) {
    double ip;
    if (ml_isnan (x)) {
        if (iptr) *iptr = x;
        return x;
    }
    if (ml_isinf (x)) {
        if (iptr) *iptr = x;
        return ml_copysign (0.0, x);
    }
    ip = ml_trunc (x);
    if (iptr) *iptr = ip;
    /* x - trunc(x) is exact for all finite doubles (Sterbenz regime). */
    if (ip == x) return ml_copysign (0.0, x);
    return x - ip;
}
ML_API double ml_scalbn (double x, int n) {
    return ml_ldexp_pure (x, n);
}
ML_API double ml_scalbln (double x, long n) {
    /* ldexp takes int; clamp huge magnitudes so the shift saturates to
     * 0 or Inf exactly as the standard requires. */
    if (n > 100000L) return ml_ldexp_pure (x, 100000);
    if (n < -100000L) return ml_ldexp_pure (x, -100000);
    return ml_ldexp_pure (x, (int) n);
}
ML_API int ml_ilogb (double x) {
    ml_fp_parts_t p;
    uint64_t s;
    int e;
    if (ml_isnan (x)) return FP_ILOGBNAN;
    if (x == 0.0 || ml_isinf (x)) return FP_ILOGB0;

    p = ml_fp_decompose (ml_fabs (x));
    s = p.sig;
    e = p.exp;
    while (s < (1ULL << 52)) {
        s <<= 1;
        e--;
    }
    return e + 52;
}
ML_API double ml_logb_fp (double x) {
    if (ml_isnan (x)) return ml_make_nan ();
    if (ml_isinf (x)) return ml_make_inf (0);
    if (x == 0.0) return ml_make_inf (1);
    return (double) ml_ilogb (x);
}
ML_API double ml_rint (double x) {
    ml_fp_parts_t p;
    int neg;
    unsigned fb;
    uint64_t ip, fmask, frac, half;

    if (ml_isnan (x) || ml_isinf (x) || x == 0.0) return x;
    neg = ml_signbit (x);
    p = ml_fp_decompose (neg ? -x : x);
    if (p.exp >= 0) return x;                      /* already integral */
    if (p.exp <= -54) return ml_copysign (0.0, x); /* |x| < 0.5 -> signed zero */

    fb = (unsigned) (-p.exp);
    ip = p.sig >> fb;
    fmask = (1ULL << fb) - 1ULL;
    frac = p.sig & fmask;
    half = 1ULL << (fb - 1);

    if (frac > half || (frac == half && (ip & 1ULL))) ip++;
    {
        double r = (double) ip;
        return neg ? -r : r;
    }
}
ML_API double ml_nearbyint (double x) {
    return ml_rint (x);
}
static long ml_to_long_checked (double x) {
    long double v = (long double) x;
    long double upper = (long double) LONG_MAX + 1.0L;
    if (!ml_isfinite (x) || v < (long double) LONG_MIN || v >= upper)
        return LONG_MIN;
    return (long) x;
}
static long long ml_to_llong_checked (double x) {
    long double v = (long double) x;
    long double upper = (long double) LLONG_MAX + 1.0L;
    if (!ml_isfinite (x) || v < (long double) LLONG_MIN || v >= upper)
        return LLONG_MIN;
    return (long long) x;
}
ML_API long ml_lrint (double x) {
    return ml_to_long_checked (ml_rint (x));
}
ML_API long long ml_llrint (double x) {
    return ml_to_llong_checked (ml_rint (x));
}
ML_API long ml_lround (double x) {
    return ml_to_long_checked (ml_round (x));
}
ML_API long long ml_llround (double x) {
    return ml_to_llong_checked (ml_round (x));
}
static int ml_u64_bit_length (uint64_t x) {
    int n = 0;
    while (x != 0u) {
        n++;
        x >>= 1;
    }
    return n;
}
/* Return seven low quotient bits, including the sign. The quotient can be
 * thousands of bits wide for valid binary64 inputs, so never convert it to a
 * native integer. Decomposition gives |x/y| = sx/sy * 2^delta exactly. */
static int ml_remquo_low7 (double x, double y) {
    ml_fp_parts_t px = ml_fp_decompose (ml_fabs (x));
    ml_fp_parts_t py = ml_fp_decompose (ml_fabs (y));
    uint64_t sx = px.sig, sy = py.sig, rem = 0u;
    int ex = px.exp, ey = py.exp, delta, nbits, i;
    int negative = ml_signbit (x) != ml_signbit (y);
    unsigned qmod = 0u;

    while (sx < (1ULL << 52)) {
        sx <<= 1;
        ex--;
    }
    while (sy < (1ULL << 52)) {
        sy <<= 1;
        ey--;
    }
    delta = ex - ey;

    if (delta >= 0) {
        nbits = ml_u64_bit_length (sx) + delta;
        for (i = nbits - 1; i >= 0; i--) {
            int src_bit = i - delta;
            unsigned in_bit = src_bit >= 0 ? (unsigned) ((sx >> src_bit) & 1u) : 0u;
            unsigned qbit;
            rem = (rem << 1) | (uint64_t) in_bit;
            qbit = rem >= sy;
            if (qbit) rem -= sy;
            qmod = ((qmod << 1) | qbit) & 0x7fu;
        }
        if ((rem << 1) > sy || ((rem << 1) == sy && (qmod & 1u)))
            qmod = (qmod + 1u) & 0x7fu;
    } else if (delta == -1 && sx > sy) {
        /* The ratio is sx/(2 sy), in [1/4,1); nearest integer is one only
         * when it is strictly above 1/2. Exact ties round to even zero. */
        qmod = 1u;
    }

    return negative ? -(int) qmod : (int) qmod;
}
ML_API double ml_remquo (double x, double y, int *quo) {
    double r;

    if (quo) *quo = 0;
    if (ml_isnan (x) || ml_isnan (y) || ml_isinf (x) || y == 0.0) return ml_make_nan ();
    if (ml_isinf (y)) return x;
    if (x == 0.0) return x;

    r = ml_remainder (x, y);
    if (quo) *quo = ml_remquo_low7 (x, y);
    return r;
}
ML_API int ml_fpclassify (double x) {
    ml_fp_parts_t p = ml_fp_decompose (x);
    switch (p.kind) {
    case ML_FP_NAN:
        return ML_CLASS_NAN;
    case ML_FP_INF:
        return ML_CLASS_INFINITE;
    case ML_FP_NORMAL:
        return ML_CLASS_NORMAL;
    case ML_FP_SUBNORMAL:
        return ML_CLASS_SUBNORMAL;
    case ML_FP_ZERO:
    default:
        return ML_CLASS_ZERO;
    }
}
ML_API int ml_isnormal (double x) {
    return ml_fpclassify (x) == ML_CLASS_NORMAL;
}
ML_API double ml_fma (double a, double b, double c) {
    return ML_FMA (a, b, c);
}
ML_API void ml_sincos (double x, double *s, double *c) {
    if (s) *s = ml_sin (x);
    if (c) *c = ml_cos (x);
}
ML_API double ml_nan (const char *tag) {
    (void) tag;
    return ml_make_nan ();
}
ML_API double ml_inf (void) {
    return ml_make_inf (0);
}
ML_API double ml_tgamma (double x) {
    return ml_gamma_new (x);
}
ML_API double ml_lgamma_r (double x, int *signgamp) {
    double v = ml_lgamma (x);
    int sg = 1;
    if (ml_isnan (x) || (ml_isinf (x) && x < 0.0)) {
        sg = 0;
    } else if (x < 0.0) {
        if (x == ml_round (x)) {
            /* Gamma has a pole here, so its sign is not a finite-value
             * property. Report zero instead of inventing one by parity. */
            sg = 0;
        } else {
            sg = ml_sinpi (x) < 0.0 ? -1 : 1;
        }
    } else if (x == 0.0 && ml_signbit (x)) {
        sg = -1;
    }
    if (signgamp) *signgamp = sg;
    return v;
}
