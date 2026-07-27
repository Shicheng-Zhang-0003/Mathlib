#include "ml_compiler.h"
#include "ml_core.h"
#include "internal/ieee_exact.h"

ML_API double ml_sqrt(double x) {
    if (x < 0.0) return ml_make_nan();
    if (x == 0.0) return x;

#if defined(__x86_64__) || defined(__i386__)
    double res;
    __asm__ ("sqrtsd %1, %0" : "=x" (res) : "x" (x));
    return res;
#else
    return __builtin_sqrt(x);
#endif
}

ML_API double ml_ldexp_pure(double x, int exp) {
    /* MATHLIB_REDP2_P0_1_LDEXP_SUBNORMAL_OVERFLOW */
    ml_fp_parts_t p = ml_fp_decompose(x);

    if (p.kind == ML_FP_ZERO || p.kind == ML_FP_INF || p.kind == ML_FP_NAN) {
        return x;
    }

    long long e = (long long)p.exp + (long long)exp;
    uint64_t sig = p.sig;

    if (sig == 0) {
        return ml_copysign(0.0, x);
    }

    /*
     * Normalize subnormal significands BEFORE applying overflow bounds.
     *
     * The old code applied the normal-number overflow bound (e > 971)
     * directly to subnormal significands, causing false overflow for
     * small sig + large positive exponent.
     *
     * Example:
     *   smallest subnormal = 2^-1074
     *   ldexp(smallest subnormal, 2046) = 2^972, which is finite.
     *
     * The old path saw:
     *   new_exp = -1074 + 2046 = 972
     *   972 > 971 -> infinity
     *
     * That was incorrect.
     */
    while (sig < (1ULL << 52) && e > -2098) {
        sig <<= 1;
        e--;
    }

    if (e > 971) {
        return ml_make_inf(p.sign);
    }

    if (e < -2098) {
        uint64_t bits = p.sign ? 0x8000000000000000ULL : 0ULL;
        double d;
        memcpy(&d, &bits, sizeof(double));
        return d;
    }

    return ml_fp_compose(sig, (int)e, p.sign);
}

ML_API double ml_frexp_pure(double x, int *exp) {
    if (!exp) return ml_make_nan();

    ml_fp_parts_t p = ml_fp_decompose(x);

    if (p.kind == ML_FP_ZERO) {
        *exp = 0;
        return x;
    }

    if (p.kind == ML_FP_INF || p.kind == ML_FP_NAN) {
        *exp = 0;
        return x;
    }

    uint64_t sig = p.sig;
    int e = p.exp;

    while (sig < (1ULL << 52)) {
        sig <<= 1;
        e--;
    }

    double m = (double)sig / 9007199254740992.0; /* 2^53 */
    *exp = e + 53;

    return p.sign ? -m : m;
}

ML_API double ml_fmod(double x, double y) {
    /* v11S AUDIT IP-1: exact fmod for all finite nonzero inputs */
    if (ml_isnan(x) || ml_isnan(y) || ml_isinf(x)) return ml_make_nan();
    if (ml_isinf(y)) return x;
    if (y == 0.0) return ml_make_nan();

    double ax = ml_fabs(x);
    double ay = ml_fabs(y);

    if (ax < ay) return x;
    if (ax == ay) return ml_copysign(0.0, x);

    ml_fp_parts_t px = ml_fp_decompose(ax);
    ml_fp_parts_t py = ml_fp_decompose(ay);

    if (px.sig == 0) return ml_copysign(0.0, x);
    if (py.sig == 0) return ml_make_nan();

    int d = px.exp - py.exp;

    if (d < 0) {
        int sd = -d;

        /*
         * For |x| >= |y| this should not happen except in exotic subnormal
         * alignments. If it does, align y to x's exponent exactly when safe.
         */
        if (sd >= 64 || py.sig > (UINT64_MAX >> sd)) {
            return ml_make_nan();
        }

        uint64_t ysig = py.sig << sd;
        uint64_t rem = px.sig % ysig;

        return ml_fp_compose(rem, px.exp, ml_signbit(x));
    }

    uint64_t rem = px.sig % py.sig;

    if (rem == 0) {
        return ml_copysign(0.0, x);
    }

    for (int i = 0; i < d; i++) {
        rem <<= 1;
        if (rem >= py.sig) rem -= py.sig;
        if (rem == 0) break;
    }

    return ml_fp_compose(rem, py.exp, ml_signbit(x));
}

ML_API double ml_round(double x) {
    /* MATHLIB_REDP2_P0_2_ROUND_CORRECT */
    if (ml_isnan(x) || ml_isinf(x)) return x;
    if (x == 0.0) return x;

    int neg = ml_signbit(x);
    ml_fp_parts_t p = ml_fp_decompose(neg ? -x : x);

    if (p.kind == ML_FP_ZERO) {
        return x;
    }

    /*
     * If exponent >= 0, the value is already an integer.
     */
    if (p.exp >= 0) {
        return x;
    }

    /*
     * If |x| < 0.5, round to signed zero.
     */
    if (p.exp <= -54) {
        return neg ? -0.0 : 0.0;
    }

    unsigned frac_bits = (unsigned)(-p.exp);

    uint64_t int_part = p.sig >> frac_bits;
    uint64_t frac_mask = (1ULL << frac_bits) - 1ULL;
    uint64_t frac_part = p.sig & frac_mask;
    uint64_t half = 1ULL << (frac_bits - 1);

    /*
     * Round half away from zero.
     *
     * This avoids the broken x +/- 0.5 trick, which can misround
     * values just below halves due to floating-point addition rounding.
     */
    if (frac_part >= half) {
        int_part++;
    }

    double r = (double)int_part;
    return neg ? -r : r;
}

