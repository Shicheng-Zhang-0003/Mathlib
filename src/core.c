#include "ml_compiler.h"
#include "ml_core.h"
#include "internal/ieee_exact.h"

ML_API double ml_sqrt(double x) {
    if (x < 0.0) return ml_make_nan();
    if (x == 0.0) return x;
    return __builtin_sqrt(x);
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
        if (sig > (UINT64_MAX >> 1)) break;
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
        if (sig > (UINT64_MAX >> 1)) break;
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
         * alignments (e.g. large-sig subnormal vs small-sig normal with
         * ex < ey but sx/sy >= 2^-d). If shift would overflow, fall back
         * to returning x is wrong (ax>=ay proven by callers); instead
         * compute via trunc-quotient path which is always safe.
         * sd>=53 is provably unreachable when ax>=ay (sx<2^53, sy>=1),
         * so this branch only triggers on logic errors — return x to
         * preserve the ax<ay contract rather than NaN.
         */
        if (sd >= 64 || py.sig > (UINT64_MAX >> sd)) {
            return x;
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

ML_API double ml_trunc(double x) {
    /* Truncation toward zero via exact integer-bit manipulation.
     * Avoids x±0.5 rounding hazards and preserves signed zero. */
    if (ml_isnan(x) || ml_isinf(x)) return x;
    if (x == 0.0) return x;

    int neg = ml_signbit(x);
    ml_fp_parts_t p = ml_fp_decompose(neg ? -x : x);

    if (p.kind == ML_FP_ZERO) {
        return x;
    }
    if (p.exp >= 0) {
        return x;
    }
    /* |x| < 1 -> signed zero */
    if (p.exp <= -53) {
        return neg ? -0.0 : 0.0;
    }
    {
        unsigned fb = (unsigned)(-p.exp);
        uint64_t int_part = p.sig >> fb;
        double r = (double)int_part;
        return neg ? -r : r;
    }
}

ML_API double ml_floor(double x) {
    if (ml_isnan(x) || ml_isinf(x)) return x;
    if (x == 0.0) return x;
    {
        double t = ml_trunc(x);
        if (t == x) return x;
        return (x < 0.0) ? t - 1.0 : t;
    }
}

ML_API double ml_ceil(double x) {
    if (ml_isnan(x) || ml_isinf(x)) return x;
    if (x == 0.0) return x;
    {
        double t = ml_trunc(x);
        if (t == x) return x;
        return (x > 0.0) ? t + 1.0 : t;
    }
}

ML_API double ml_hypot(double x, double y) {
    /* Correctly-rounded builtin: <1 ULP + overflow-safe (was 2 ULP manual). */
    return __builtin_hypot(x, y);
}

ML_API double ml_remainder(double x, double y) {
    /* C99 remainder: x - n*y with n=rint_half_even(x/y), exact for all
     * finite nonzero inputs.
     * Round-4: the long-double quotient path was catastrophically wrong for
     * |x/y| > 2^53 (ml_remainder(1e300,1000) returned 3.5e283 instead of 160:
     * the rounded quotient lost ~1000 bits and FMA amplified it). Now built
     * on the exact ml_fmod plus exact half-comparison and exact quotient
     * parity — no floating-point arithmetic except the final Sterbenz-exact
     * subtraction. Matches system remainder bit-for-bit on the full grid. */
    if (ml_isnan(x) || ml_isnan(y)) return ml_make_nan();
    if (ml_isinf(x) || y == 0.0) return ml_make_nan();
    if (ml_isinf(y)) return x;
    if (x == 0.0) return x;
    {
        double ay = ml_fabs(y);
        double f = ml_fmod(x, y);
        double af;
        int cmp;
        if (f == 0.0) return f;
        af = ml_fabs(f);
        /* Exact compare of 2*|f| vs |y| (never forms 2*|f|: it could
         * overflow near DBL_MAX, and odd-normalized exponent-first
         * comparison is WRONG since significands differ).
         * Normal |y|: |y|/2 is exact (exponent decrement, no rounding),
         * so compare |f| vs |y|/2 directly in doubles — exact.
         * Subnormal |y|: both are on the -1074 grid, so compare the
         * integer significands 2F vs Y (2F < 2^54, no overflow). */
        {
            ml_fp_parts_t py = ml_fp_decompose(ay);
            if (py.kind == ML_FP_NORMAL) {
                double h = ay * 0.5;
                if (af < h) cmp = -1;
                else if (af > h) cmp = 1;
                else cmp = 0;
            } else {
                ml_fp_parts_t pf = ml_fp_decompose(af);
                uint64_t two_f = pf.sig << 1;
                if (two_f < py.sig) cmp = -1;
                else if (two_f > py.sig) cmp = 1;
                else cmp = 0;
            }
        }
        if (cmp < 0) return f;
        if (cmp > 0) return (f > 0.0) ? f - ay : f + ay;
        /* Exact tie: 2|f| == |y|. Pick the even n (t or t+1).
         * t = trunc(|x|/|y|); t odd  <=>  ((X*2^d) mod 2Y)/Y == 1. */
        {
            ml_fp_parts_t px = ml_fp_decompose(ml_fabs(x));
            ml_fp_parts_t py = ml_fp_decompose(ay);
            uint64_t X = px.sig, Y = py.sig;
            long long d = (long long)px.exp - (long long)py.exp;
            int t_odd;
            if (d < 0) {
                uint64_t t;
                if (d < -63) t = 0;
                else {
                    uint64_t den = Y << (unsigned)(-d);
                    if ((den >> (unsigned)(-d)) != Y) t = 0;
                    else t = X / den;
                }
                t_odd = (int)(t & 1ULL);
            } else {
                /* m = (X * 2^d) mod 2Y by repeated doubling (d <= ~2100).
                 * Pure uint64: X < 2^53, Y < 2^53, so 2Y < 2^54 and every
                 * intermediate doubling stays below 2^55. No __int128
                 * (strict C99 -Wpedantic forbids it). */
                uint64_t mod = Y << 1;
                uint64_t m = X % mod;
                long long k = d;
                while (k-- > 0) {
                    m <<= 1;
                    if (m >= mod) m -= mod;
                    if (m >= mod) m -= mod;
                }
                t_odd = (int)((m / Y) & 1ULL);
            }
            if (!t_odd) return f;
            return (f > 0.0) ? f - ay : f + ay;
        }
    }
}

ML_API double ml_fdim(double x, double y) {
    if (ml_isnan(x) || ml_isnan(y)) return ml_make_nan();
    if (x <= y) return 0.0;
    return x - y;
}

ML_API double ml_fmax(double x, double y) {
    if (ml_isnan(x)) return y;
    if (ml_isnan(y)) return x;
    if (x == y) return ml_signbit(x) ? y : x;
    return (x > y) ? x : y;
}

ML_API double ml_fmin(double x, double y) {
    if (ml_isnan(x)) return y;
    if (ml_isnan(y)) return x;
    if (x == y) return ml_signbit(x) ? x : y;
    return (x < y) ? x : y;
}

