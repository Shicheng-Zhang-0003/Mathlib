#include "ml_compiler.h"
#include "ml_bigint.h"
#include <limits.h>
#include <string.h>

#define BI_LIMBS ML_BIGINT_LIMBS
/* The public struct is stack-allocated by callers, so validate its bounded
 * representation before indexing any limb. */
static int bi_valid (const ml_bigint_t *a) {
    int i;
    if (a == NULL || a->n < 0 || a->n > BI_LIMBS ||
        (a->sign != 0 && a->sign != 1)) return 0;
    if (a->n == 0) {
        if (a->sign != 0) return 0;
    } else if (a->limb [a->n - 1] == 0u) {
        return 0;
    }
    for (i = a->n; i < BI_LIMBS; i++)
        if (a->limb [i] != 0u) return 0;
    return 1;
}
static void bi_trim (ml_bigint_t *a) {
    while (a->n > 0 && a->limb [a->n - 1] == 0u) a->n--;
    if (a->n == 0) a->sign = 0;
}
static int bi_cmp_mag (const ml_bigint_t *a, const ml_bigint_t *b) {
    int i;
    if (a->n != b->n) return a->n < b->n ? -1 : 1;
    for (i = a->n - 1; i >= 0; i--) {
        if (a->limb [i] != b->limb [i]) return a->limb [i] < b->limb [i] ? -1 : 1;
    }
    return 0;
}
static ml_status_t bi_add_mag (ml_bigint_t *r, const ml_bigint_t *a, const ml_bigint_t *b) {
    int i, m = a->n > b->n ? a->n : b->n;
    uint32_t carry = 0u;
    for (i = 0; i < m; i++) {
        uint32_t av = i < a->n ? a->limb [i] : 0u;
        uint32_t bv = i < b->n ? b->limb [i] : 0u;
        uint64_t t = (uint64_t) av + (uint64_t) bv + (uint64_t) carry;
        r->limb [i] = (uint32_t) t;
        carry = (uint32_t) (t >> 32);
    }
    r->n = m;
    if (carry) {
        if (r->n >= BI_LIMBS) return ML_ERR_OVERFLOW;
        r->limb [r->n++] = carry;
    }
    r->sign = 0;
    bi_trim (r);
    return ML_SUCCESS;
}
/* r = |a| - |b|, requires |a| >= |b| */
static ml_status_t bi_sub_mag (ml_bigint_t *r, const ml_bigint_t *a, const ml_bigint_t *b) {
    int i;
    uint32_t borrow = 0u;
    for (i = 0; i < a->n; i++) {
        uint32_t av = a->limb [i];
        uint32_t bv = i < b->n ? b->limb [i] : 0u;
        uint64_t t = (uint64_t) av - (uint64_t) bv - (uint64_t) borrow;
        r->limb [i] = (uint32_t) t;
        borrow = (t >> 32) ? 1u : 0u;
    }
    r->n = a->n;
    r->sign = 0;
    bi_trim (r);
    return ML_SUCCESS;
}
static ml_status_t bi_mul_small (ml_bigint_t *a, uint32_t m) {
    int i;
    uint64_t carry = 0;
    for (i = 0; i < a->n; i++) {
        uint64_t t = (uint64_t) a->limb [i] * (uint64_t) m + carry;
        a->limb [i] = (uint32_t) t;
        carry = t >> 32;
    }
    while (carry) {
        if (a->n >= BI_LIMBS) return ML_ERR_OVERFLOW;
        a->limb [a->n++] = (uint32_t) carry;
        carry >>= 32;
    }
    bi_trim (a);
    return ML_SUCCESS;
}
static ml_status_t bi_add_small (ml_bigint_t *a, uint32_t v) {
    int i = 0;
    uint64_t carry = v;
    while (carry) {
        uint64_t t;
        if (i >= a->n) {
            if (i >= BI_LIMBS) return ML_ERR_OVERFLOW;
            a->limb [i] = 0u;
            a->n = i + 1;
        }
        t = (uint64_t) a->limb [i] + carry;
        a->limb [i] = (uint32_t) t;
        carry = t >> 32;
        i++;
    }
    bi_trim (a);
    return ML_SUCCESS;
}
static uint32_t bi_div_small (ml_bigint_t *a, uint32_t d) {
    int i;
    uint64_t rem = 0;
    for (i = a->n - 1; i >= 0; i--) {
        uint64_t cur = (rem << 32) | (uint64_t) a->limb [i];
        a->limb [i] = (uint32_t) (cur / (uint64_t) d);
        rem = cur % (uint64_t) d;
    }
    bi_trim (a);
    return (uint32_t) rem;
}
static int bi_shl1 (ml_bigint_t *a) {
    int i;
    uint32_t carry = 0u;
    for (i = 0; i < a->n; i++) {
        uint32_t v = a->limb [i];
        a->limb [i] = (uint32_t) ((v << 1) | carry);
        carry = v >> 31;
    }
    if (carry) {
        if (a->n >= BI_LIMBS) return 0;
        a->limb [a->n++] = carry;
    }
    return 1;
}
static void bi_set_bit (ml_bigint_t *q, int i) {
    int l = i >> 5, b = i & 31;
    if (l >= BI_LIMBS) return;
    q->limb [l] |= (uint32_t) (1u << b);
    if (l >= q->n) { q->n = l + 1; }
}
static int bitlen32 (uint32_t v) {
    int c = 0;
    while (v) {
        c++;
        v >>= 1;
    }
    return c;
}
ML_API void ml_bigint_init (ml_bigint_t *a) {
    if (a == NULL) return;
    memset (a->limb, 0, sizeof (a->limb));
    a->n = 0;
    a->sign = 0;
}
ML_API ml_status_t ml_bigint_set_u64 (ml_bigint_t *a, uint64_t v) {
    if (a == NULL) return ML_ERR_INVALID_ARG;
    ml_bigint_init (a);
    if (v == 0) return ML_SUCCESS;
    a->limb [0] = (uint32_t) v;
    a->n = 1;
    if ((v >> 32) != 0) {
        a->limb [1] = (uint32_t) (v >> 32);
        a->n = 2;
    }
    return ML_SUCCESS;
}
ML_API ml_status_t ml_bigint_set_i64 (ml_bigint_t *a, int64_t v) {
    uint64_t mag;
    if (a == NULL) return ML_ERR_INVALID_ARG;
    if (v < 0) {
        mag = (uint64_t) (-(v + 1)) + 1u;
        (void) ml_bigint_set_u64 (a, mag);
        a->sign = a->n ? 1 : 0;
    } else {
        (void) ml_bigint_set_u64 (a, (uint64_t) v);
    }
    return ML_SUCCESS;
}
ML_API ml_status_t ml_bigint_copy (ml_bigint_t *r, const ml_bigint_t *a) {
    if (r == NULL || !bi_valid (a)) return ML_ERR_INVALID_ARG;
    *r = *a;
    return ML_SUCCESS;
}
static int digit_val (char c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'z') return c - 'a' + 10;
    if (c >= 'A' && c <= 'Z') return c - 'A' + 10;
    return -1;
}
static int bi_is_space (char c) {
    return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\v' || c == '\f';
}
ML_API ml_status_t ml_bigint_from_str (ml_bigint_t *a, const char *s, int base) {
    ml_bigint_t parsed;
    size_t i = 0u;
    int neg = 0, saw_digit = 0, trailing_space = 0;
    if (a == NULL || s == NULL || base < 2 || base > 36) return ML_ERR_INVALID_ARG;
    ml_bigint_init (&parsed);
    while (bi_is_space (s [i])) i++;
    if (s [i] == '+' || s [i] == '-') {
        neg = (s [i] == '-');
        i++;
    }
    for (; s [i]; i++) {
        int d = digit_val (s [i]);
        ml_status_t st;
        if (bi_is_space (s [i])) {
            trailing_space = 1;
            continue;
        }
        if (trailing_space || d < 0 || d >= base) return ML_ERR_INVALID_ARG;
        st = bi_mul_small (&parsed, (uint32_t) base);
        if (st != ML_SUCCESS) return st;
        st = bi_add_small (&parsed, (uint32_t) d);
        if (st != ML_SUCCESS) return st;
        saw_digit = 1;
    }
    if (!saw_digit) return ML_ERR_INVALID_ARG;
    if (neg && parsed.n) parsed.sign = 1;
    *a = parsed;
    return ML_SUCCESS;
}
static const char *bi_digits = "0123456789abcdefghijklmnopqrstuvwxyz";
ML_API int ml_bigint_to_str (const ml_bigint_t *a, char *buf, size_t bufsz, int base) {
    ml_bigint_t tmp;
    char rev [4200];
    int len = 0, i;
    size_t need;
    if (!bi_valid (a) || !buf || base < 2 || base > 36) return -1;
    if (a->n == 0) {
        if (bufsz < 2) return -1;
        buf [0] = '0';
        buf [1] = '\0';
        return 1;
    }
    tmp = *a;
    tmp.sign = 0;
    while (tmp.n > 0) {
        uint32_t d = bi_div_small (&tmp, (uint32_t) base);
        if (len >= (int) sizeof (rev)) return -1;
        rev [len++] = bi_digits [d];
    }
    need = (size_t) len + (a->sign ? 1u : 0u) + 1u;
    if (need > bufsz) return -1;
    i = 0;
    if (a->sign) buf [i++] = '-';
    while (len > 0) buf [i++] = rev [--len];
    buf [i] = '\0';
    return i;
}
ML_API int ml_bigint_cmp (const ml_bigint_t *a, const ml_bigint_t *b) {
    if (!bi_valid (a) || !bi_valid (b)) return INT_MIN;
    if (a->n == 0 && b->n == 0) return 0;
    if (a->sign != b->sign) return a->sign ? -1 : 1;
    {
        int c = bi_cmp_mag (a, b);
        return a->sign ? -c : c;
    }
}
ML_API int ml_bigint_is_zero (const ml_bigint_t *a) { return bi_valid (a) && a->n == 0; }
ML_API int ml_bigint_bitlen (const ml_bigint_t *a) {
    if (!bi_valid (a)) return -1;
    if (a->n == 0) return 0;
    return (a->n - 1) * 32 + bitlen32 (a->limb [a->n - 1]);
}
ML_API double ml_bigint_to_double (const ml_bigint_t *a) {
    double v = 0.0;
    int i;
    if (!bi_valid (a)) return ml_make_nan ();
    for (i = a->n - 1; i >= 0; i--) v = v * 4294967296.0 + (double) a->limb [i];
    return a->sign ? -v : v;
}
ML_API ml_status_t ml_bigint_abs (ml_bigint_t *r, const ml_bigint_t *a) {
    if (r == NULL || !bi_valid (a)) return ML_ERR_INVALID_ARG;
    *r = *a;
    r->sign = 0;
    return ML_SUCCESS;
}
ML_API ml_status_t ml_bigint_neg (ml_bigint_t *r, const ml_bigint_t *a) {
    ml_bigint_t t;
    if (r == NULL || !bi_valid (a)) return ML_ERR_INVALID_ARG;
    t = *a;
    if (t.n) t.sign = t.sign ? 0 : 1;
    *r = t;
    return ML_SUCCESS;
}
ML_API ml_status_t ml_bigint_add (ml_bigint_t *r, const ml_bigint_t *a, const ml_bigint_t *b) {
    ml_bigint_t t;
    ml_status_t st;
    if (r == NULL || !bi_valid (a) || !bi_valid (b)) return ML_ERR_INVALID_ARG;
    ml_bigint_init (&t);
    if (a->sign == b->sign) {
        st = bi_add_mag (&t, a, b);
        if (st != ML_SUCCESS) return st;
        t.sign = a->n ? a->sign : 0;
    } else {
        int c = bi_cmp_mag (a, b);
        if (c == 0) {
            ml_bigint_init (r);
            return ML_SUCCESS;
        }
        if (c > 0) {
            (void) bi_sub_mag (&t, a, b);
            t.sign = a->sign;
        } else {
            (void) bi_sub_mag (&t, b, a);
            t.sign = b->sign;
        }
    }
    *r = t;
    return ML_SUCCESS;
}
ML_API ml_status_t ml_bigint_sub (ml_bigint_t *r, const ml_bigint_t *a, const ml_bigint_t *b) {
    if (r == NULL || !bi_valid (a) || !bi_valid (b)) return ML_ERR_INVALID_ARG;
    ml_bigint_t nb = *b;
    if (nb.n) nb.sign = nb.sign ? 0 : 1;
    return ml_bigint_add (r, a, &nb);
}
ML_API ml_status_t ml_bigint_mul (ml_bigint_t *r, const ml_bigint_t *a, const ml_bigint_t *b) {
    ml_bigint_t t;
    int i, j;
    if (r == NULL || !bi_valid (a) || !bi_valid (b)) return ML_ERR_INVALID_ARG;
    if (a->n == 0 || b->n == 0) {
        ml_bigint_init (r);
        return ML_SUCCESS;
    }
    /* A product of n- and m-limb values needs at most n+m limbs, and at
     * least n+m-1.  The n+m==capacity+1 case can still fit, so let the
     * carry propagation decide whether it really overflows. */
    if (a->n + b->n > BI_LIMBS + 1) return ML_ERR_OVERFLOW;
    ml_bigint_init (&t);
    for (i = 0; i < a->n; i++) {
        uint64_t carry = 0;
        uint32_t ai = a->limb [i];
        for (j = 0; j < b->n; j++) {
            uint64_t prod = (uint64_t) ai * (uint64_t) b->limb [j];
            uint64_t cur = (uint64_t) t.limb [i + j] + (prod & 0xFFFFFFFFULL) + carry;
            t.limb [i + j] = (uint32_t) cur;
            carry = (prod >> 32) + (cur >> 32);
        }
        {
            int k = i + b->n;
            while (carry) {
                uint64_t cur;
                if (k >= BI_LIMBS) return ML_ERR_OVERFLOW;
                cur = (uint64_t) t.limb [k] + carry;
                t.limb [k] = (uint32_t) cur;
                carry = cur >> 32;
                k++;
            }
            t.n = k > t.n ? k : t.n;
        }
    }
    t.n = a->n + b->n;
    if (t.n > BI_LIMBS) t.n = BI_LIMBS;
    bi_trim (&t);
    t.sign = (a->sign != b->sign) && t.n ? 1 : 0;
    *r = t;
    return ML_SUCCESS;
}
ML_API ml_status_t ml_bigint_div (ml_bigint_t *q, ml_bigint_t *rem,
                                  const ml_bigint_t *a, const ml_bigint_t *b) {
    ml_bigint_t qq, rr, am, bm;
    int i, bl;
    if (!bi_valid (a) || !bi_valid (b) || (q == NULL && rem == NULL) ||
        (q != NULL && q == rem) || b->n == 0) return ML_ERR_INVALID_ARG;
    am = *a;
    am.sign = 0;
    bm = *b;
    bm.sign = 0;
    ml_bigint_init (&qq);
    ml_bigint_init (&rr);
    if (bi_cmp_mag (&am, &bm) < 0) {
        rr = am;
    } else {
        bl = ml_bigint_bitlen (&am);
        for (i = bl - 1; i >= 0; i--) {
            if (!bi_shl1 (&rr)) return ML_ERR_OVERFLOW;
            if ((am.limb [i >> 5] >> (i & 31)) & 1u) {
                rr.limb [0] |= 1u;
                if (rr.n == 0) rr.n = 1;
            }
            if (bi_cmp_mag (&rr, &bm) >= 0) {
                (void) bi_sub_mag (&rr, &rr, &bm);
                bi_set_bit (&qq, i);
            }
        }
        bi_trim (&qq);
    }
    bi_trim (&rr);
    qq.sign = (a->sign != b->sign) && qq.n ? 1 : 0;
    rr.sign = a->sign && rr.n ? 1 : 0;
    if (q) *q = qq;
    if (rem) *rem = rr;
    return ML_SUCCESS;
}
ML_API ml_status_t ml_bigint_mod (ml_bigint_t *r, const ml_bigint_t *a, const ml_bigint_t *b) {
    ml_bigint_t rr, modulus;
    if (r == NULL || !bi_valid (a) || !bi_valid (b) || b->n == 0)
        return ML_ERR_INVALID_ARG;
    ml_status_t st = ml_bigint_div (NULL, &rr, a, b);
    if (st != ML_SUCCESS) return st;
    modulus = *b;
    modulus.sign = 0;
    /* make non-negative */
    if (rr.sign) {
        st = ml_bigint_add (&rr, &rr, &modulus);
        if (st != ML_SUCCESS) return st;
    }
    *r = rr;
    return ML_SUCCESS;
}
ML_API ml_status_t ml_bigint_shl (ml_bigint_t *r, const ml_bigint_t *a, unsigned bits) {
    ml_bigint_t t;
    unsigned limb_shift = bits >> 5, bit_shift = bits & 31u;
    uint64_t out_bits;
    int i;
    if (r == NULL || !bi_valid (a)) return ML_ERR_INVALID_ARG;
    if (a->n == 0) {
        ml_bigint_init (r);
        return ML_SUCCESS;
    }
    if (bits == 0u) {
        *r = *a;
        return ML_SUCCESS;
    }
    out_bits = (uint64_t) (unsigned) ml_bigint_bitlen (a) + (uint64_t) bits;
    if (out_bits > (uint64_t) BI_LIMBS * 32u) return ML_ERR_OVERFLOW;
    ml_bigint_init (&t);
    for (i = 0; i < a->n; i++) {
        uint64_t v = (uint64_t) a->limb [i] << bit_shift;
        int k = i + (int) limb_shift;
        t.limb [k] |= (uint32_t) (v & 0xFFFFFFFFULL);
        if (k + 1 < BI_LIMBS) t.limb [k + 1] |= (uint32_t) (v >> 32);
    }
    t.n = a->n + (int) limb_shift + 1;
    if (t.n > BI_LIMBS) t.n = BI_LIMBS;
    t.sign = a->sign;
    bi_trim (&t);
    *r = t;
    return ML_SUCCESS;
}
ML_API ml_status_t ml_bigint_shr (ml_bigint_t *r, const ml_bigint_t *a, unsigned bits) {
    ml_bigint_t t;
    unsigned limb_shift = bits >> 5, bit_shift = bits & 31u;
    int i;
    if (r == NULL || !bi_valid (a)) return ML_ERR_INVALID_ARG;
    ml_bigint_init (&t);
    if ((int) limb_shift >= a->n) {
        *r = t;
        return ML_SUCCESS;
    }
    for (i = 0; i + (int) limb_shift < a->n; i++) {
        int src = i + (int) limb_shift;
        uint64_t v = (uint64_t) a->limb [src] >> bit_shift;
        if (bit_shift && src + 1 < a->n) {
            v |= (uint64_t) a->limb [src + 1] << (32 - bit_shift);
        }
        t.limb [i] = (uint32_t) v;
    }
    t.n = a->n - (int) limb_shift;
    t.sign = a->sign;
    bi_trim (&t);
    *r = t;
    return ML_SUCCESS;
}
ML_API ml_status_t ml_bigint_pow (ml_bigint_t *r, const ml_bigint_t *base, uint64_t exp) {
    ml_bigint_t result, b, tmp;
    ml_status_t st;
    if (r == NULL || !bi_valid (base)) return ML_ERR_INVALID_ARG;
    int odd = (int) (exp & 1u);
    (void) ml_bigint_set_u64 (&result, 1);
    b = *base;
    b.sign = 0;
    while (exp) {
        if (exp & 1u) {
            st = ml_bigint_mul (&tmp, &result, &b);
            if (st != ML_SUCCESS) return st;
            result = tmp;
        }
        exp >>= 1;
        if (exp) {
            st = ml_bigint_mul (&tmp, &b, &b);
            if (st != ML_SUCCESS) return st;
            b = tmp;
        }
    }
    result.sign = (base->sign && odd && result.n) ? 1 : 0;
    *r = result;
    return ML_SUCCESS;
}
/* Add two reduced non-negative values without forming a possibly 4097-bit
 * intermediate.  Preconditions: a,b < mod and mod > 0. */
static ml_status_t bi_add_mod (ml_bigint_t *r, const ml_bigint_t *a,
                               const ml_bigint_t *b, const ml_bigint_t *mod) {
    ml_bigint_t threshold, t;
    ml_status_t st = ml_bigint_sub (&threshold, mod, b);
    if (st != ML_SUCCESS) return st;
    ml_bigint_init (&t);
    if (bi_cmp_mag (a, &threshold) >= 0) {
        (void) bi_sub_mag (&t, a, &threshold);
    } else {
        st = bi_add_mag (&t, a, b);
        if (st != ML_SUCCESS) return st;
    }
    t.sign = 0;
    *r = t;
    return ML_SUCCESS;
}
/* Binary modular multiplication keeps every intermediate below mod, so it
 * works even when the full product would exceed the fixed limb capacity. */
static ml_status_t bi_mul_mod (ml_bigint_t *r, const ml_bigint_t *a,
                               const ml_bigint_t *b, const ml_bigint_t *mod) {
    ml_bigint_t acc, x = *a, y = *b, t;
    ml_status_t st;
    ml_bigint_init (&acc);
    while (y.n) {
        if (y.limb [0] & 1u) {
            st = bi_add_mod (&t, &acc, &x, mod);
            if (st != ML_SUCCESS) return st;
            acc = t;
        }
        st = ml_bigint_shr (&y, &y, 1u);
        if (st != ML_SUCCESS) return st;
        if (y.n) {
            st = bi_add_mod (&t, &x, &x, mod);
            if (st != ML_SUCCESS) return st;
            x = t;
        }
    }
    *r = acc;
    return ML_SUCCESS;
}
ML_API ml_status_t ml_bigint_powmod (ml_bigint_t *r, const ml_bigint_t *base,
                                     const ml_bigint_t *exp, const ml_bigint_t *mod) {
    ml_bigint_t result, b, tmp;
    int bl, i;
    ml_status_t st;
    if (r == NULL || !bi_valid (base) || !bi_valid (exp) || !bi_valid (mod) ||
        exp->sign || mod->sign || mod->n == 0) return ML_ERR_INVALID_ARG;
    (void) ml_bigint_set_u64 (&tmp, 1);
    st = ml_bigint_mod (&result, &tmp, mod);
    if (st != ML_SUCCESS) return st;
    st = ml_bigint_mod (&b, base, mod);
    if (st != ML_SUCCESS) return st;
    bl = ml_bigint_bitlen (exp);
    for (i = 0; i < bl; i++) {
        if ((exp->limb [i >> 5] >> (i & 31)) & 1u) {
            st = bi_mul_mod (&tmp, &result, &b, mod);
            if (st != ML_SUCCESS) return st;
            result = tmp;
        }
        st = bi_mul_mod (&tmp, &b, &b, mod);
        if (st != ML_SUCCESS) return st;
        b = tmp;
    }
    *r = result;
    return ML_SUCCESS;
}
ML_API ml_status_t ml_bigint_gcd (ml_bigint_t *r, const ml_bigint_t *a, const ml_bigint_t *b) {
    ml_bigint_t x, y, t;
    if (r == NULL || !bi_valid (a) || !bi_valid (b)) return ML_ERR_INVALID_ARG;
    x = *a;
    y = *b;
    x.sign = 0;
    y.sign = 0;
    while (y.n) {
        ml_status_t st = ml_bigint_mod (&t, &x, &y);
        if (st != ML_SUCCESS) return st;
        x = y;
        y = t;
    }
    *r = x;
    return ML_SUCCESS;
}
ML_API ml_status_t ml_bigint_lcm (ml_bigint_t *r, const ml_bigint_t *a, const ml_bigint_t *b) {
    ml_bigint_t g, aa, bb, q;
    ml_status_t st;
    if (r == NULL || !bi_valid (a) || !bi_valid (b)) return ML_ERR_INVALID_ARG;
    if (a->n == 0 || b->n == 0) {
        ml_bigint_init (r);
        return ML_SUCCESS;
    }
    st = ml_bigint_gcd (&g, a, b);
    if (st != ML_SUCCESS) return st;
    aa = *a;
    aa.sign = 0;
    bb = *b;
    bb.sign = 0;
    st = ml_bigint_div (&q, NULL, &aa, &g);
    if (st != ML_SUCCESS) return st;
    st = ml_bigint_mul (&q, &q, &bb);
    if (st != ML_SUCCESS) return st;
    q.sign = 0;
    *r = q;
    return ML_SUCCESS;
}
ML_API ml_status_t ml_bigint_modinv (ml_bigint_t *r, const ml_bigint_t *a, const ml_bigint_t *m) {
    ml_bigint_t old_r, rr, old_s, s, q, tmp, prod;
    ml_status_t st;
    if (r == NULL || !bi_valid (a) || !bi_valid (m) || m->sign || m->n == 0)
        return ML_ERR_INVALID_ARG;
    st = ml_bigint_mod (&old_r, a, m);
    if (st != ML_SUCCESS) return st;
    rr = *m;
    (void) ml_bigint_set_u64 (&old_s, 1);
    ml_bigint_init (&s);
    old_r.sign = 0;
    rr.sign = 0;
    while (rr.n) {
        st = ml_bigint_div (&q, &tmp, &old_r, &rr);
        if (st != ML_SUCCESS) return st;
        old_r = rr;
        rr = tmp;
        st = ml_bigint_mul (&prod, &q, &s);
        if (st != ML_SUCCESS) return st;
        st = ml_bigint_sub (&tmp, &old_s, &prod);
        if (st != ML_SUCCESS) return st;
        old_s = s;
        s = tmp;
    }
    if (!(old_r.n == 1 && old_r.limb [0] == 1u)) return ML_ERR_UNSUPPORTED;
    st = ml_bigint_mod (r, &old_s, m);
    if (st != ML_SUCCESS) return st;
    return ML_SUCCESS;
}
ML_API ml_status_t ml_bigint_factorial (ml_bigint_t *r, unsigned n) {
    ml_bigint_t a;
    unsigned i;
    if (r == NULL) return ML_ERR_INVALID_ARG;
    (void) ml_bigint_set_u64 (&a, 1);
    for (i = 2; i <= n; i++) {
        ml_status_t st = bi_mul_small (&a, i);
        if (st != ML_SUCCESS) return st;
    }
    *r = a;
    return ML_SUCCESS;
}
ML_API ml_status_t ml_bigint_fib (ml_bigint_t *r, unsigned n) {
    ml_bigint_t a, b, t;
    unsigned i;
    if (r == NULL) return ML_ERR_INVALID_ARG;
    ml_bigint_init (&a);
    (void) ml_bigint_set_u64 (&b, 1);
    if (n == 0) {
        ml_bigint_init (r);
        return ML_SUCCESS;
    }
    for (i = 1; i < n; i++) {
        ml_status_t st = ml_bigint_add (&t, &a, &b);
        if (st != ML_SUCCESS) return st;
        a = b;
        b = t;
    }
    *r = b;
    return ML_SUCCESS;
}
static uint64_t bi_xorshift (uint64_t *state) {
    uint64_t x = *state;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    *state = x;
    return x * 2685821657736338717ULL;
}
ML_API ml_status_t ml_bigint_rand (ml_bigint_t *r, unsigned bits, uint64_t *state) {
    unsigned limbs;
    unsigned i;
    if (r == NULL || state == NULL || bits == 0u || bits > BI_LIMBS * 32u)
        return ML_ERR_INVALID_ARG;
    limbs = (bits + 31u) / 32u;
    /* xorshift64* has a zero absorbing state; give seed zero a defined,
     * reproducible nonzero stream. This generator is not cryptographic. */
    if (*state == 0u) *state = 0x9E3779B97F4A7C15ULL;
    ml_bigint_init (r);
    for (i = 0; i < limbs; i++) r->limb [i] = (uint32_t) bi_xorshift (state);
    r->n = (int) limbs;
    {
        unsigned top = bits & 31u;
        if (top) r->limb [limbs - 1] &= (uint32_t) ((1u << top) - 1u);
    }
    bi_trim (r);
    return ML_SUCCESS;
}
ML_API int ml_bigint_is_probable_prime (const ml_bigint_t *a, int rounds, uint64_t *state) {
    static const uint32_t small [] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    ml_bigint_t n, one, two, three, nm1, d, x, base, tmp, span;
    int i, s = 0, j;
    unsigned sm;
    if (!bi_valid (a) || a->sign || a->n <= 0 || rounds <= 0 || rounds > 256)
        return 0;
    if (a->n == 1 && a->limb [0] < 2u) return 0;
    if (a->n == 1 && (a->limb [0] == 2u || a->limb [0] == 3u)) return 1;
    if ((a->limb [0] & 1u) == 0u) return 0;

    /* Exact trial division for 32-bit inputs. Testing only the first twelve
     * primes is insufficient: 41^2, for example, would otherwise pass. */
    if (a->n == 1) {
        uint32_t v = a->limb [0];
        for (sm = 0; sm < sizeof (small) / sizeof (small [0]); sm++) {
            if (v == small [sm]) return 1;
            if (v % small [sm] == 0u) return 0;
        }
        for (uint32_t d32 = 41u; d32 <= v / d32; d32 += 2u)
            if (v % d32 == 0u) return 0;
        return 1;
    }
    for (sm = 0; sm < sizeof (small) / sizeof (small [0]); sm++) {
        if (a->n == 1 && a->limb [0] == small [sm]) return 1;
        {
            ml_bigint_t rem;
            (void) ml_bigint_set_u64 (&tmp, small [sm]);
            if (ml_bigint_div (NULL, &rem, a, &tmp) != ML_SUCCESS) return 0;
            if (rem.n == 0) return 0;
        }
    }
    if (state == NULL) return 0;
    n = *a;
    (void) ml_bigint_set_u64 (&one, 1);
    (void) ml_bigint_set_u64 (&two, 2);
    (void) ml_bigint_set_u64 (&three, 3);
    (void) ml_bigint_sub (&nm1, &n, &one);
    d = nm1;
    while (d.n && (d.limb [0] & 1u) == 0u) {
        (void) ml_bigint_shr (&d, &d, 1);
        s++;
    }
    (void) ml_bigint_sub (&span, &n, &three);
    if (span.n == 0) return 0;
    for (j = 0; j < rounds; j++) {
        ml_status_t st;
        int tries;
        for (tries = 0; tries < 128; tries++) {
            st = ml_bigint_rand (&base, (unsigned) ml_bigint_bitlen (&span), state);
            if (st != ML_SUCCESS) return 0;
            if (bi_cmp_mag (&base, &span) < 0) break;
        }
        if (tries == 128) return 0;
        (void) ml_bigint_add (&base, &base, &two);
        st = ml_bigint_powmod (&x, &base, &d, &n);
        if (st != ML_SUCCESS) return 0;
        if (x.n == 0) return 0;
        if (x.n == 1 && x.limb [0] == 1u) continue;
        if (ml_bigint_cmp (&x, &nm1) == 0) continue;
        {
            int composite = 1;
            for (i = 1; i < s; i++) {
                st = ml_bigint_mul (&tmp, &x, &x);
                if (st != ML_SUCCESS) return 0;
                st = ml_bigint_mod (&x, &tmp, &n);
                if (st != ML_SUCCESS) return 0;
                if (ml_bigint_cmp (&x, &nm1) == 0) {
                    composite = 0;
                    break;
                }
                if (x.n == 1 && x.limb [0] == 1u) break;
            }
            if (composite) return 0;
        }
    }
    return 1;
}
