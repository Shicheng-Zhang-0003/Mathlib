#include "ml_exact.h"
#include "ml_core.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
static void mpz_zero (ml_mpz_t *a) {
    a->size = 0u;
    a->sign = 0;
}
static void mpz_normalize (ml_mpz_t *a) {
    while (a->size > 0u && a->limb [a->size - 1u] == 0u) a->size--;
    if (a->size == 0u) a->sign = 0;
}
static void mpz_swap (ml_mpz_t *a, ml_mpz_t *b) {
    ml_mpz_t t = *a;
    *a = *b;
    *b = t;
}
static ml_status_t mpz_reserve (ml_mpz_t *a, size_t needed) {
    size_t cap;
    uint32_t *p;
    if (needed <= a->capacity) return ML_SUCCESS;
    if (needed > SIZE_MAX / sizeof (uint32_t)) return ML_ERR_OVERFLOW;
    cap = a->capacity ? a->capacity : 4u;
    while (cap < needed) {
        if (cap > SIZE_MAX / 2u) {
            cap = needed;
            break;
        }
        cap *= 2u;
    }
    if (cap > SIZE_MAX / sizeof (uint32_t)) cap = needed;
    p = (uint32_t *) realloc (a->limb, cap * sizeof (uint32_t));
    if (p == NULL) return ML_ERR_WORKSPACE;
    a->limb = p;
    a->capacity = cap;
    return ML_SUCCESS;
}
static void mpz_commit (ml_mpz_t *out, ml_mpz_t *tmp) {
    ml_mpz_clear (out);
    *out = *tmp;
    ml_mpz_init (tmp);
}
void ml_mpz_init (ml_mpz_t *a) {
    if (a == NULL) return;
    a->limb = NULL;
    a->size = 0u;
    a->capacity = 0u;
    a->sign = 0;
}
void ml_mpz_clear (ml_mpz_t *a) {
    if (a == NULL) return;
    free (a->limb);
    ml_mpz_init (a);
}
int ml_mpz_is_valid (const ml_mpz_t *a) {
    if (a == NULL || (a->sign != 0 && a->sign != 1) || a->size > a->capacity)
        return 0;
    if ((a->capacity == 0u) != (a->limb == NULL)) return 0;
    if (a->size == 0u) return a->sign == 0;
    return a->limb [a->size - 1u] != 0u;
}
int ml_mpz_is_zero (const ml_mpz_t *a) {
    return ml_mpz_is_valid (a) && a->size == 0u;
}
int ml_mpz_sign (const ml_mpz_t *a) {
    if (!ml_mpz_is_valid (a) || a->size == 0u) return 0;
    return a->sign ? -1 : 1;
}
size_t ml_mpz_bit_length (const ml_mpz_t *a) {
    uint32_t top;
    size_t bits = 0u;
    if (!ml_mpz_is_valid (a) || a->size == 0u) return 0u;
    top = a->limb [a->size - 1u];
    while (top != 0u) {
        bits++;
        top >>= 1;
    }
    if (a->size - 1u > (SIZE_MAX - bits) / 32u) return SIZE_MAX;
    return 32u * (a->size - 1u) + bits;
}
static int mpz_cmp_mag (const ml_mpz_t *a, const ml_mpz_t *b) {
    size_t i;
    if (a->size != b->size) return a->size < b->size ? -1 : 1;
    i = a->size;
    while (i > 0u) {
        i--;
        if (a->limb [i] != b->limb [i]) return a->limb [i] < b->limb [i] ? -1 : 1;
    }
    return 0;
}
int ml_mpz_compare (const ml_mpz_t *a, const ml_mpz_t *b) {
    int c;
    if (!ml_mpz_is_valid (a) || !ml_mpz_is_valid (b)) return INT_MIN;
    if (a->size == 0u && b->size == 0u) return 0;
    if (a->sign != b->sign) return a->sign ? -1 : 1;
    c = mpz_cmp_mag (a, b);
    return a->sign ? -c : c;
}
static ml_status_t mpz_copy_raw (ml_mpz_t *r, const ml_mpz_t *a) {
    ml_status_t st;
    if (r == a) return ML_SUCCESS;
    st = mpz_reserve (r, a->size);
    if (st != ML_SUCCESS) return st;
    if (a->size) memcpy (r->limb, a->limb, a->size * sizeof (uint32_t));
    r->size = a->size;
    r->sign = a->sign;
    return ML_SUCCESS;
}
ml_status_t ml_mpz_copy (ml_mpz_t *r, const ml_mpz_t *a) {
    if (!ml_mpz_is_valid (r) || !ml_mpz_is_valid (a)) return ML_ERR_INVALID_ARG;
    return mpz_copy_raw (r, a);
}
static ml_status_t mpz_set_u64_raw (ml_mpz_t *r, uint64_t v) {
    ml_status_t st = mpz_reserve (r, v > UINT32_MAX ? 2u : (v ? 1u : 0u));
    if (st != ML_SUCCESS) return st;
    mpz_zero (r);
    if (v) {
        r->limb [0] = (uint32_t) v;
        r->size = 1u;
        if (v >> 32) {
            r->limb [1] = (uint32_t) (v >> 32);
            r->size = 2u;
        }
    }
    return ML_SUCCESS;
}
ml_status_t ml_mpz_set_u64 (ml_mpz_t *r, uint64_t value) {
    ml_mpz_t t;
    ml_status_t st;
    if (!ml_mpz_is_valid (r)) return ML_ERR_INVALID_ARG;
    ml_mpz_init (&t);
    st = mpz_set_u64_raw (&t, value);
    if (st == ML_SUCCESS) mpz_commit (r, &t);
    ml_mpz_clear (&t);
    return st;
}
ml_status_t ml_mpz_set_i64 (ml_mpz_t *r, int64_t value) {
    uint64_t mag;
    ml_mpz_t t;
    ml_status_t st;
    if (!ml_mpz_is_valid (r)) return ML_ERR_INVALID_ARG;
    mag = value < 0 ? (uint64_t) (-(value + 1)) + 1u : (uint64_t) value;
    ml_mpz_init (&t);
    st = mpz_set_u64_raw (&t, mag);
    if (st == ML_SUCCESS) {
        t.sign = value < 0 && t.size != 0u;
        mpz_commit (r, &t);
    }
    ml_mpz_clear (&t);
    return st;
}
static ml_status_t mpz_mul_small (ml_mpz_t *a, uint32_t m) {
    size_t i;
    uint64_t carry = 0u;
    ml_status_t st;
    if (a->size == 0u || m == 0u) {
        mpz_zero (a);
        return ML_SUCCESS;
    }
    st = mpz_reserve (a, a->size + 1u);
    if (st != ML_SUCCESS) return st;
    for (i = 0u; i < a->size; i++) {
        uint64_t z = (uint64_t) a->limb [i] * m + carry;
        a->limb [i] = (uint32_t) z;
        carry = z >> 32;
    }
    if (carry) a->limb [a->size++] = (uint32_t) carry;
    return ML_SUCCESS;
}
static ml_status_t mpz_add_small (ml_mpz_t *a, uint32_t v) {
    size_t i = 0u;
    uint64_t carry = v;
    ml_status_t st;
    if (v == 0u) return ML_SUCCESS;
    if (a->sign) return ML_ERR_INTERNAL;
    st = mpz_reserve (a, a->size + 1u);
    if (st != ML_SUCCESS) return st;
    while (carry != 0u) {
        uint64_t z;
        if (i == a->size) a->limb [a->size++] = 0u;
        z = (uint64_t) a->limb [i] + carry;
        a->limb [i] = (uint32_t) z;
        carry = z >> 32;
        i++;
    }
    return ML_SUCCESS;
}
static int digit_value (char c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'z') return c - 'a' + 10;
    if (c >= 'A' && c <= 'Z') return c - 'A' + 10;
    return -1;
}
ml_status_t ml_mpz_set_string (ml_mpz_t *r, const char *text, unsigned base) {
    ml_mpz_t t;
    size_t i = 0u;
    int negative = 0, any = 0;
    ml_status_t st = ML_SUCCESS;
    if (!ml_mpz_is_valid (r) || text == NULL || base < 2u || base > 36u)
        return ML_ERR_INVALID_ARG;
    if (text [i] == '+' || text [i] == '-') {
        negative = text [i] == '-';
        i++;
    }
    ml_mpz_init (&t);
    for (; text [i] != '\0'; i++) {
        int d = digit_value (text [i]);
        if (d < 0 || (unsigned) d >= base) {
            st = ML_ERR_INVALID_ARG;
            goto done;
        }
        st = mpz_mul_small (&t, (uint32_t) base);
        if (st != ML_SUCCESS) goto done;
        st = mpz_add_small (&t, (uint32_t) d);
        if (st != ML_SUCCESS) goto done;
        any = 1;
    }
    if (!any) {
        st = ML_ERR_INVALID_ARG;
        goto done;
    }
    t.sign = negative && t.size != 0u;
    mpz_commit (r, &t);
done:
    ml_mpz_clear (&t);
    return st;
}
static uint32_t mpz_div_small_inplace (ml_mpz_t *a, uint32_t d) {
    size_t i = a->size;
    uint64_t rem = 0u;
    while (i > 0u) {
        uint64_t z;
        i--;
        z = (rem << 32) | a->limb [i];
        a->limb [i] = (uint32_t) (z / d);
        rem = z % d;
    }
    mpz_normalize (a);
    return (uint32_t) rem;
}
int ml_mpz_to_string (const ml_mpz_t *a, char *buffer, size_t size, unsigned base) {
    static const char digits [] = "0123456789abcdefghijklmnopqrstuvwxyz";
    ml_mpz_t t;
    char *rev;
    size_t cap, n = 0u, need, i;
    if (!ml_mpz_is_valid (a) || buffer == NULL || base < 2u || base > 36u)
        return -1;
    if (a->size == 0u) {
        if (size < 2u) return -1;
        buffer [0] = '0';
        buffer [1] = '\0';
        return 1;
    }
    cap = ml_mpz_bit_length (a);
    if (cap > SIZE_MAX - 3u) return -1;
    cap += 3u;
    rev = (char *) malloc (cap);
    if (rev == NULL) return -1;
    ml_mpz_init (&t);
    if (mpz_copy_raw (&t, a) != ML_SUCCESS) {
        free (rev);
        return -1;
    }
    t.sign = 0;
    while (t.size != 0u) {
        if (n >= cap) {
            ml_mpz_clear (&t);
            free (rev);
            return -1;
        }
        rev [n++] = digits [mpz_div_small_inplace (&t, (uint32_t) base)];
    }
    need = n + (a->sign ? 1u : 0u) + 1u;
    if (need > size || n > (size_t) INT_MAX - (a->sign ? 1u : 0u)) {
        ml_mpz_clear (&t);
        free (rev);
        return -1;
    }
    i = 0u;
    if (a->sign) buffer [i++] = '-';
    while (n > 0u) buffer [i++] = rev [--n];
    buffer [i] = '\0';
    ml_mpz_clear (&t);
    free (rev);
    return (int) i;
}
static ml_status_t mpz_add_mag (ml_mpz_t *r, const ml_mpz_t *a,
                                const ml_mpz_t *b) {
    size_t i, n = a->size > b->size ? a->size : b->size;
    uint64_t carry = 0u;
    ml_status_t st;
    if (n == SIZE_MAX) return ML_ERR_OVERFLOW;
    st = mpz_reserve (r, n + 1u);
    if (st != ML_SUCCESS) return st;
    for (i = 0u; i < n; i++) {
        uint64_t av = i < a->size ? a->limb [i] : 0u;
        uint64_t bv = i < b->size ? b->limb [i] : 0u;
        uint64_t z = av + bv + carry;
        r->limb [i] = (uint32_t) z;
        carry = z >> 32;
    }
    r->size = n;
    if (carry) r->limb [r->size++] = (uint32_t) carry;
    r->sign = 0;
    return ML_SUCCESS;
}
/* Requires |a| >= |b|. */
static ml_status_t mpz_sub_mag (ml_mpz_t *r, const ml_mpz_t *a,
                                const ml_mpz_t *b) {
    size_t i;
    uint64_t borrow = 0u;
    ml_status_t st = mpz_reserve (r, a->size);
    if (st != ML_SUCCESS) return st;
    for (i = 0u; i < a->size; i++) {
        uint64_t av = a->limb [i];
        uint64_t bv = (i < b->size ? b->limb [i] : 0u) + borrow;
        r->limb [i] = (uint32_t) (av - bv);
        borrow = av < bv;
    }
    r->size = a->size;
    r->sign = 0;
    mpz_normalize (r);
    return ML_SUCCESS;
}
ml_status_t ml_mpz_add (ml_mpz_t *r, const ml_mpz_t *a, const ml_mpz_t *b) {
    ml_mpz_t t;
    ml_status_t st;
    int c;
    if (!ml_mpz_is_valid (r) || !ml_mpz_is_valid (a) || !ml_mpz_is_valid (b))
        return ML_ERR_INVALID_ARG;
    ml_mpz_init (&t);
    if (a->sign == b->sign) {
        st = mpz_add_mag (&t, a, b);
        t.sign = t.size ? a->sign : 0;
    } else {
        c = mpz_cmp_mag (a, b);
        if (c >= 0) {
            st = mpz_sub_mag (&t, a, b);
            t.sign = t.size ? a->sign : 0;
        } else {
            st = mpz_sub_mag (&t, b, a);
            t.sign = t.size ? b->sign : 0;
        }
    }
    if (st == ML_SUCCESS) mpz_commit (r, &t);
    ml_mpz_clear (&t);
    return st;
}
ml_status_t ml_mpz_sub (ml_mpz_t *r, const ml_mpz_t *a, const ml_mpz_t *b) {
    ml_mpz_t neg;
    ml_status_t st;
    if (!ml_mpz_is_valid (r) || !ml_mpz_is_valid (a) || !ml_mpz_is_valid (b))
        return ML_ERR_INVALID_ARG;
    ml_mpz_init (&neg);
    st = mpz_copy_raw (&neg, b);
    if (st == ML_SUCCESS && neg.size) neg.sign = !neg.sign;
    if (st == ML_SUCCESS) st = ml_mpz_add (r, a, &neg);
    ml_mpz_clear (&neg);
    return st;
}
ml_status_t ml_mpz_mul (ml_mpz_t *r, const ml_mpz_t *a, const ml_mpz_t *b) {
    ml_mpz_t t;
    ml_status_t st;
    size_t i, j, n;
    if (!ml_mpz_is_valid (r) || !ml_mpz_is_valid (a) || !ml_mpz_is_valid (b))
        return ML_ERR_INVALID_ARG;
    ml_mpz_init (&t);
    if (a->size == 0u || b->size == 0u) {
        mpz_commit (r, &t);
        return ML_SUCCESS;
    }
    if (a->size > SIZE_MAX - b->size) return ML_ERR_OVERFLOW;
    n = a->size + b->size;
    st = mpz_reserve (&t, n);
    if (st != ML_SUCCESS) {
        ml_mpz_clear (&t);
        return st;
    }
    memset (t.limb, 0, n * sizeof (uint32_t));
    t.size = n;
    for (i = 0u; i < a->size; i++) {
        uint64_t carry = 0u;
        for (j = 0u; j < b->size; j++) {
            size_t k = i + j;
            uint64_t z = (uint64_t) a->limb [i] * b->limb [j] + t.limb [k] + carry;
            t.limb [k] = (uint32_t) z;
            carry = z >> 32;
        }
        {
            size_t k = i + b->size;
            while (carry != 0u && k < n) {
                uint64_t z = (uint64_t) t.limb [k] + carry;
                t.limb [k] = (uint32_t) z;
                carry = z >> 32;
                k++;
            }
        }
    }
    t.sign = a->sign != b->sign;
    mpz_normalize (&t);
    mpz_commit (r, &t);
    ml_mpz_clear (&t);
    return ML_SUCCESS;
}
static ml_status_t mpz_shift_one_inplace (ml_mpz_t *a) {
    size_t i;
    uint32_t carry = 0u;
    ml_status_t st;
    if (a->size == 0u) return ML_SUCCESS;
    if (a->size == SIZE_MAX) return ML_ERR_OVERFLOW;
    st = mpz_reserve (a, a->size + 1u);
    if (st != ML_SUCCESS) return st;
    for (i = 0u; i < a->size; i++) {
        uint32_t v = a->limb [i];
        a->limb [i] = (v << 1) | carry;
        carry = v >> 31;
    }
    if (carry) a->limb [a->size++] = carry;
    return ML_SUCCESS;
}
static int mpz_get_bit (const ml_mpz_t *a, size_t bit) {
    size_t i = bit / 32u;
    return i < a->size ? (int) ((a->limb [i] >> (bit % 32u)) & 1u) : 0;
}
static ml_status_t mpz_set_bit (ml_mpz_t *a, size_t bit) {
    size_t i = bit / 32u;
    ml_status_t st;
    if (i == SIZE_MAX) return ML_ERR_OVERFLOW;
    st = mpz_reserve (a, i + 1u);
    if (st != ML_SUCCESS) return st;
    while (a->size <= i) a->limb [a->size++] = 0u;
    a->limb [i] |= (uint32_t) (1u << (bit % 32u));
    return ML_SUCCESS;
}
static ml_status_t mpz_sub_mag_inplace (ml_mpz_t *a, const ml_mpz_t *b) {
    size_t i;
    uint64_t borrow = 0u;
    if (mpz_cmp_mag (a, b) < 0) return ML_ERR_INTERNAL;
    for (i = 0u; i < a->size; i++) {
        uint64_t av = a->limb [i];
        uint64_t bv = (i < b->size ? b->limb [i] : 0u) + borrow;
        a->limb [i] = (uint32_t) (av - bv);
        borrow = av < bv;
    }
    mpz_normalize (a);
    return ML_SUCCESS;
}
ml_status_t ml_mpz_divrem (ml_mpz_t *quotient, ml_mpz_t *remainder,
                           const ml_mpz_t *a, const ml_mpz_t *b) {
    ml_mpz_t q, rem;
    size_t bits, i;
    ml_status_t st = ML_SUCCESS;
    if (!ml_mpz_is_valid (quotient) || !ml_mpz_is_valid (remainder) ||
        quotient == remainder || !ml_mpz_is_valid (a) || !ml_mpz_is_valid (b))
        return ML_ERR_INVALID_ARG;
    if (b->size == 0u) return ML_ERR_SINGULAR;
    ml_mpz_init (&q);
    ml_mpz_init (&rem);
    bits = ml_mpz_bit_length (a);
    if (bits == SIZE_MAX) {
        st = ML_ERR_OVERFLOW;
        goto done;
    }
    if (bits != 0u) {
        size_t limbs = bits / 32u + (bits % 32u != 0u);
        st = mpz_reserve (&q, limbs);
        if (st != ML_SUCCESS) goto done;
        memset (q.limb, 0, limbs * sizeof (uint32_t));
        q.size = limbs;
    }
    i = bits;
    while (i > 0u) {
        i--;
        st = mpz_shift_one_inplace (&rem);
        if (st != ML_SUCCESS) goto done;
        if (mpz_get_bit (a, i)) {
            st = mpz_add_small (&rem, 1u);
            if (st != ML_SUCCESS) goto done;
        }
        if (mpz_cmp_mag (&rem, b) >= 0) {
            st = mpz_sub_mag_inplace (&rem, b);
            if (st != ML_SUCCESS) goto done;
            st = mpz_set_bit (&q, i);
            if (st != ML_SUCCESS) goto done;
        }
    }
    mpz_normalize (&q);
    q.sign = q.size ? (a->sign != b->sign) : 0;
    rem.sign = rem.size ? a->sign : 0;
    mpz_commit (quotient, &q);
    mpz_commit (remainder, &rem);
done:
    ml_mpz_clear (&q);
    ml_mpz_clear (&rem);
    return st;
}
ml_status_t ml_mpz_gcd (ml_mpz_t *r, const ml_mpz_t *a, const ml_mpz_t *b) {
    ml_mpz_t x, y, q, rem;
    ml_status_t st;
    if (!ml_mpz_is_valid (r) || !ml_mpz_is_valid (a) || !ml_mpz_is_valid (b))
        return ML_ERR_INVALID_ARG;
    ml_mpz_init (&x);
    ml_mpz_init (&y);
    ml_mpz_init (&q);
    ml_mpz_init (&rem);
    st = mpz_copy_raw (&x, a);
    if (st == ML_SUCCESS) st = mpz_copy_raw (&y, b);
    if (st != ML_SUCCESS) goto done;
    x.sign = 0;
    y.sign = 0;
    while (y.size != 0u) {
        st = ml_mpz_divrem (&q, &rem, &x, &y);
        if (st != ML_SUCCESS) goto done;
        mpz_swap (&x, &y);
        mpz_swap (&y, &rem);
        ml_mpz_clear (&q);
    }
    mpz_commit (r, &x);
done:
    ml_mpz_clear (&x);
    ml_mpz_clear (&y);
    ml_mpz_clear (&q);
    ml_mpz_clear (&rem);
    return st;
}
ml_status_t ml_mpz_shift_left (ml_mpz_t *r, const ml_mpz_t *a, size_t bits) {
    ml_mpz_t t;
    size_t words, rembits, n, i;
    uint32_t carry = 0u;
    ml_status_t st;
    if (!ml_mpz_is_valid (r) || !ml_mpz_is_valid (a)) return ML_ERR_INVALID_ARG;
    ml_mpz_init (&t);
    if (a->size == 0u) {
        mpz_commit (r, &t);
        return ML_SUCCESS;
    }
    words = bits / 32u;
    rembits = bits % 32u;
    if (words > SIZE_MAX - a->size - 1u) return ML_ERR_OVERFLOW;
    n = a->size + words + (rembits != 0u);
    st = mpz_reserve (&t, n);
    if (st != ML_SUCCESS) {
        ml_mpz_clear (&t);
        return st;
    }
    memset (t.limb, 0, n * sizeof (uint32_t));
    for (i = 0u; i < a->size; i++) {
        uint64_t z = ((uint64_t) a->limb [i] << rembits) | carry;
        t.limb [i + words] = (uint32_t) z;
        carry = (uint32_t) (z >> 32);
    }
    if (carry) t.limb [a->size + words] = carry;
    t.size = n;
    t.sign = a->sign;
    mpz_normalize (&t);
    mpz_commit (r, &t);
    ml_mpz_clear (&t);
    return ML_SUCCESS;
}
double ml_mpz_to_double (const ml_mpz_t *a) {
    size_t start, i;
    long double v = 0.0L;
    int exp2;
    double d;
    if (!ml_mpz_is_valid (a)) return ml_make_nan ();
    if (a->size == 0u) return 0.0;
    start = a->size > 3u ? a->size - 3u : 0u;
    i = a->size;
    while (i > start) {
        i--;
        v = v * 4294967296.0L + (long double) a->limb [i];
    }
    if (start > (size_t) (INT_MAX / 32)) return a->sign ? -ml_make_inf (0) : ml_make_inf (0);
    exp2 = (int) (start * 32u);
    v = ldexpl (v, exp2);
    d = (double) v;
    return a->sign ? -d : d;
}
/* ------------------------------ canonical rationals ------------------- */

void ml_mpq_init (ml_mpq_t *q) {
    if (q == NULL) return;
    ml_mpz_init (&q->numerator);
    ml_mpz_init (&q->denominator);
}
void ml_mpq_clear (ml_mpq_t *q) {
    if (q == NULL) return;
    ml_mpz_clear (&q->numerator);
    ml_mpz_clear (&q->denominator);
}
int ml_mpq_is_valid (const ml_mpq_t *q) {
    return q != NULL && ml_mpz_is_valid (&q->numerator) &&
           ml_mpz_is_valid (&q->denominator) && q->denominator.size != 0u &&
           q->denominator.sign == 0;
}
static int mpq_initialized (const ml_mpq_t *q) {
    return q != NULL && ml_mpz_is_valid (&q->numerator) &&
           ml_mpz_is_valid (&q->denominator);
}
static void mpq_commit (ml_mpq_t *out, ml_mpq_t *tmp) {
    ml_mpq_clear (out);
    *out = *tmp;
    ml_mpq_init (tmp);
}
static ml_status_t mpz_div_exact (ml_mpz_t *out, const ml_mpz_t *a,
                                  const ml_mpz_t *b) {
    ml_mpz_t q, rem;
    ml_status_t st;
    ml_mpz_init (&q);
    ml_mpz_init (&rem);
    st = ml_mpz_divrem (&q, &rem, a, b);
    if (st == ML_SUCCESS && rem.size != 0u) st = ML_ERR_INTERNAL;
    if (st == ML_SUCCESS) mpz_commit (out, &q);
    ml_mpz_clear (&q);
    ml_mpz_clear (&rem);
    return st;
}
static ml_status_t mpq_normalize (ml_mpq_t *q) {
    ml_mpz_t g, n, d;
    ml_status_t st;
    ml_mpz_init (&g);
    ml_mpz_init (&n);
    ml_mpz_init (&d);
    if (q->denominator.size == 0u) {
        st = ML_ERR_SINGULAR;
        goto done;
    }
    if (q->denominator.sign) {
        q->denominator.sign = 0;
        if (q->numerator.size) q->numerator.sign = !q->numerator.sign;
    }
    if (q->numerator.size == 0u) {
        st = mpz_set_u64_raw (&q->denominator, 1u);
        goto done;
    }
    st = ml_mpz_gcd (&g, &q->numerator, &q->denominator);
    if (st != ML_SUCCESS) goto done;
    st = mpz_div_exact (&n, &q->numerator, &g);
    if (st != ML_SUCCESS) goto done;
    st = mpz_div_exact (&d, &q->denominator, &g);
    if (st != ML_SUCCESS) goto done;
    mpz_swap (&q->numerator, &n);
    mpz_swap (&q->denominator, &d);
done:
    ml_mpz_clear (&g);
    ml_mpz_clear (&n);
    ml_mpz_clear (&d);
    return st;
}
static ml_status_t mpz_set_pow10 (ml_mpz_t *out, uint32_t exponent) {
    ml_mpz_t accumulator, base, temporary;
    ml_status_t st;
    ml_mpz_init (&accumulator);
    ml_mpz_init (&base);
    ml_mpz_init (&temporary);
    st = ml_mpz_set_u64 (&accumulator, 1u);
    if (st == ML_SUCCESS) st = ml_mpz_set_u64 (&base, 10u);
    while (st == ML_SUCCESS && exponent != 0u) {
        if ((exponent & 1u) != 0u) {
            st = ml_mpz_mul (&temporary, &accumulator, &base);
            if (st == ML_SUCCESS) mpz_swap (&accumulator, &temporary);
        }
        exponent >>= 1;
        if (st == ML_SUCCESS && exponent != 0u) {
            st = ml_mpz_mul (&temporary, &base, &base);
            if (st == ML_SUCCESS) mpz_swap (&base, &temporary);
        }
    }
    if (st == ML_SUCCESS) mpz_commit (out, &accumulator);
    ml_mpz_clear (&accumulator);
    ml_mpz_clear (&base);
    ml_mpz_clear (&temporary);
    return st;
}
ml_status_t ml_mpq_set_i64 (ml_mpq_t *r, int64_t value) {
    ml_mpq_t t;
    ml_status_t st;
    if (r == NULL || !ml_mpz_is_valid (&r->numerator) ||
        !ml_mpz_is_valid (&r->denominator)) return ML_ERR_INVALID_ARG;
    ml_mpq_init (&t);
    st = ml_mpz_set_i64 (&t.numerator, value);
    if (st == ML_SUCCESS) st = ml_mpz_set_u64 (&t.denominator, 1u);
    if (st == ML_SUCCESS) mpq_commit (r, &t);
    ml_mpq_clear (&t);
    return st;
}
ml_status_t ml_mpq_set_ratio (ml_mpq_t *r, const ml_mpz_t *numerator,
                              const ml_mpz_t *denominator) {
    ml_mpq_t t;
    ml_status_t st;
    if (r == NULL || !ml_mpz_is_valid (&r->numerator) ||
        !ml_mpz_is_valid (&r->denominator) ||
        !ml_mpz_is_valid (numerator) || !ml_mpz_is_valid (denominator))
        return ML_ERR_INVALID_ARG;
    ml_mpq_init (&t);
    st = mpz_copy_raw (&t.numerator, numerator);
    if (st == ML_SUCCESS) st = mpz_copy_raw (&t.denominator, denominator);
    if (st == ML_SUCCESS) st = mpq_normalize (&t);
    if (st == ML_SUCCESS) mpq_commit (r, &t);
    ml_mpq_clear (&t);
    return st;
}
ml_status_t ml_mpq_set_decimal (ml_mpq_t *r, const char *text) {
    enum { MAX_SIGNIFICAND_DIGITS = 10000,
           MAX_EXPONENT = 10000,
           MAX_TEXT_LENGTH = 20004 };
    ml_mpz_t numerator, denominator, power, scaled;
    ml_mpz_t *numerator_to_use = &numerator;
    size_t length = 0u, i = 0u, significand_digits = 0u, fractional_digits = 0u;
    uint32_t exponent_magnitude = 0u;
    int negative = 0, exponent_negative = 0, saw_point = 0;
    int saw_exponent_digit = 0;
    int64_t exponent, decimal_scale;
    ml_status_t st = ML_SUCCESS;
    if (!mpq_initialized (r) || text == NULL) return ML_ERR_INVALID_ARG;
    while (length <= (size_t) MAX_TEXT_LENGTH && text [length] != '\0') length++;
    if (length == 0u) return ML_ERR_INVALID_ARG;
    if (length > (size_t) MAX_TEXT_LENGTH) return ML_ERR_UNSUPPORTED;

    ml_mpz_init (&numerator);
    ml_mpz_init (&denominator);
    ml_mpz_init (&power);
    ml_mpz_init (&scaled);
    if (text [i] == '+' || text [i] == '-') {
        negative = text [i] == '-';
        i++;
    }
    while (i < length && text [i] != 'e' && text [i] != 'E') {
        unsigned digit;
        if (text [i] >= '0' && text [i] <= '9') {
            digit = (unsigned) (text [i] - '0');
            if (significand_digits == (size_t) MAX_SIGNIFICAND_DIGITS) {
                st = ML_ERR_UNSUPPORTED;
                goto done;
            }
            st = mpz_mul_small (&numerator, 10u);
            if (st == ML_SUCCESS) st = mpz_add_small (&numerator, (uint32_t) digit);
            if (st != ML_SUCCESS) goto done;
            significand_digits++;
            if (saw_point) fractional_digits++;
        } else if (text [i] == '.' && !saw_point) {
            saw_point = 1;
        } else {
            st = ML_ERR_INVALID_ARG;
            goto done;
        }
        i++;
    }
    if (significand_digits == 0u) {
        st = ML_ERR_INVALID_ARG;
        goto done;
    }
    if (i < length) {
        uint32_t limit = (uint32_t) MAX_EXPONENT;
        i++;
        if (i < length && (text [i] == '+' || text [i] == '-')) {
            exponent_negative = text [i] == '-';
            i++;
        }
        if (i == length) {
            st = ML_ERR_INVALID_ARG;
            goto done;
        }
        for (; i < length; i++) {
            uint32_t digit;
            if (text [i] < '0' || text [i] > '9') {
                st = ML_ERR_INVALID_ARG;
                goto done;
            }
            digit = (uint32_t) (text [i] - '0');
            if (exponent_magnitude > (limit - digit) / 10u) {
                st = ML_ERR_UNSUPPORTED;
                goto done;
            }
            exponent_magnitude = exponent_magnitude * 10u + digit;
            saw_exponent_digit = 1;
        }
        if (!saw_exponent_digit) {
            st = ML_ERR_INVALID_ARG;
            goto done;
        }
    }
    if (numerator.size == 0u) {
        st = ml_mpq_set_i64 (r, 0);
        goto done;
    }
    numerator.sign = negative;
    exponent = exponent_negative ? -(int64_t) exponent_magnitude : (int64_t) exponent_magnitude;
    decimal_scale = (int64_t) fractional_digits - exponent;
    if (decimal_scale > 20000 || decimal_scale < -20000) {
        st = ML_ERR_UNSUPPORTED;
        goto done;
    }
    if (decimal_scale > 0) {
        st = mpz_set_pow10 (&denominator, (uint32_t) decimal_scale);
    } else {
        st = ml_mpz_set_u64 (&denominator, 1u);
        if (st == ML_SUCCESS && decimal_scale < 0) {
            st = mpz_set_pow10 (&power, (uint32_t) (-decimal_scale));
            if (st == ML_SUCCESS) st = ml_mpz_mul (&scaled, &numerator, &power);
            if (st == ML_SUCCESS) numerator_to_use = &scaled;
        }
    }
    if (st == ML_SUCCESS) st = ml_mpq_set_ratio (r, numerator_to_use, &denominator);
done:
    ml_mpz_clear (&numerator);
    ml_mpz_clear (&denominator);
    ml_mpz_clear (&power);
    ml_mpz_clear (&scaled);
    return st;
}
ml_status_t ml_mpq_copy (ml_mpq_t *r, const ml_mpq_t *q) {
    ml_mpq_t t;
    ml_status_t st;
    if (!mpq_initialized (r) || !ml_mpq_is_valid (q))
        return ML_ERR_INVALID_ARG;
    if (r == q) return ML_SUCCESS;
    ml_mpq_init (&t);
    st = mpz_copy_raw (&t.numerator, &q->numerator);
    if (st == ML_SUCCESS) st = mpz_copy_raw (&t.denominator, &q->denominator);
    if (st == ML_SUCCESS) mpq_commit (r, &t);
    ml_mpq_clear (&t);
    return st;
}
static ml_status_t mpq_binary (ml_mpq_t *r, const ml_mpq_t *a,
                               const ml_mpq_t *b, int op) {
    ml_mpz_t x, y, n, d;
    ml_mpq_t t;
    ml_status_t st = ML_SUCCESS;
    if (!mpq_initialized (r) || !ml_mpq_is_valid (a) ||
        !ml_mpq_is_valid (b)) return ML_ERR_INVALID_ARG;
    ml_mpz_init (&x);
    ml_mpz_init (&y);
    ml_mpz_init (&n);
    ml_mpz_init (&d);
    ml_mpq_init (&t);
    if (op == 3 && b->numerator.size == 0u) {
        st = ML_ERR_SINGULAR;
        goto done;
    }
    if (op == 0 || op == 1) {
        st = ml_mpz_mul (&x, &a->numerator, &b->denominator);
        if (st != ML_SUCCESS) goto done;
        st = ml_mpz_mul (&y, &b->numerator, &a->denominator);
        if (st != ML_SUCCESS) goto done;
        st = op == 0 ? ml_mpz_add (&n, &x, &y) : ml_mpz_sub (&n, &x, &y);
        if (st != ML_SUCCESS) goto done;
        st = ml_mpz_mul (&d, &a->denominator, &b->denominator);
    } else if (op == 2) {
        st = ml_mpz_mul (&n, &a->numerator, &b->numerator);
        if (st == ML_SUCCESS) st = ml_mpz_mul (&d, &a->denominator, &b->denominator);
    } else {
        st = ml_mpz_mul (&n, &a->numerator, &b->denominator);
        if (st == ML_SUCCESS) st = ml_mpz_mul (&d, &a->denominator, &b->numerator);
    }
    if (st != ML_SUCCESS) goto done;
    st = mpz_copy_raw (&t.numerator, &n);
    if (st == ML_SUCCESS) st = mpz_copy_raw (&t.denominator, &d);
    if (st == ML_SUCCESS) st = mpq_normalize (&t);
    if (st == ML_SUCCESS) mpq_commit (r, &t);
done:
    ml_mpz_clear (&x);
    ml_mpz_clear (&y);
    ml_mpz_clear (&n);
    ml_mpz_clear (&d);
    ml_mpq_clear (&t);
    return st;
}
ml_status_t ml_mpq_add (ml_mpq_t *r, const ml_mpq_t *a, const ml_mpq_t *b) {
    return mpq_binary (r, a, b, 0);
}
ml_status_t ml_mpq_sub (ml_mpq_t *r, const ml_mpq_t *a, const ml_mpq_t *b) {
    return mpq_binary (r, a, b, 1);
}
ml_status_t ml_mpq_mul (ml_mpq_t *r, const ml_mpq_t *a, const ml_mpq_t *b) {
    return mpq_binary (r, a, b, 2);
}
ml_status_t ml_mpq_div (ml_mpq_t *r, const ml_mpq_t *a, const ml_mpq_t *b) {
    return mpq_binary (r, a, b, 3);
}
ml_status_t ml_mpq_compare (const ml_mpq_t *a, const ml_mpq_t *b,
                            int *comparison) {
    ml_mpz_t x, y;
    ml_status_t st;
    if (!ml_mpq_is_valid (a) || !ml_mpq_is_valid (b) || comparison == NULL)
        return ML_ERR_INVALID_ARG;
    ml_mpz_init (&x);
    ml_mpz_init (&y);
    st = ml_mpz_mul (&x, &a->numerator, &b->denominator);
    if (st == ML_SUCCESS) st = ml_mpz_mul (&y, &b->numerator, &a->denominator);
    if (st == ML_SUCCESS) *comparison = ml_mpz_compare (&x, &y);
    ml_mpz_clear (&x);
    ml_mpz_clear (&y);
    return st;
}
static ml_status_t mpq_from_double (ml_mpq_t *q, double x) {
    ml_mpz_t n, d, shifted;
    int e;
    double f;
    uint64_t mant;
    ml_status_t st;
    ml_mpz_init (&n);
    ml_mpz_init (&d);
    ml_mpz_init (&shifted);
    if (!ml_isfinite (x)) {
        st = ML_ERR_INVALID_ARG;
        goto done;
    }
    if (x == 0.0) {
        st = ml_mpq_set_i64 (q, 0);
        goto done;
    }
    f = frexp (ml_fabs (x), &e);
    mant = (uint64_t) ldexp (f, 53);
    st = ml_mpz_set_u64 (&n, mant);
    if (st != ML_SUCCESS) goto done;
    if (x < 0.0) n.sign = 1;
    st = ml_mpz_set_u64 (&d, 1u);
    if (st != ML_SUCCESS) goto done;
    if (e >= 53) {
        st = ml_mpz_shift_left (&shifted, &n, (size_t) (e - 53));
        if (st != ML_SUCCESS) goto done;
        st = mpz_copy_raw (&n, &shifted);
    } else {
        st = ml_mpz_shift_left (&shifted, &d, (size_t) (53 - e));
        if (st != ML_SUCCESS) goto done;
        st = mpz_copy_raw (&d, &shifted);
    }
    if (st == ML_SUCCESS) st = ml_mpq_set_ratio (q, &n, &d);
done:
    ml_mpz_clear (&n);
    ml_mpz_clear (&d);
    ml_mpz_clear (&shifted);
    return st;
}
static uint64_t double_order_key (double x) {
    uint64_t bits;
    memcpy (&bits, &x, sizeof (bits));
    return (bits & UINT64_C (0x8000000000000000))
               ? ~bits
               : (bits ^ UINT64_C (0x8000000000000000));
}
static double double_from_order_key (uint64_t key) {
    uint64_t bits = (key & UINT64_C (0x8000000000000000))
                        ? (key ^ UINT64_C (0x8000000000000000))
                        : ~key;
    double x;
    memcpy (&x, &bits, sizeof (x));
    return x;
}
static ml_status_t mpq_compare_double (const ml_mpq_t *q, double x, int *cmp) {
    ml_mpq_t d;
    ml_status_t st;
    ml_mpq_init (&d);
    st = mpq_from_double (&d, x);
    if (st == ML_SUCCESS) st = ml_mpq_compare (q, &d, cmp);
    ml_mpq_clear (&d);
    return st;
}
ml_status_t ml_mpq_to_interval (const ml_mpq_t *q, ml_interval_t *out) {
    uint64_t lo, hi;
    double d;
    int cmp;
    ml_status_t st;
    ml_interval_t r;
    if (!ml_mpq_is_valid (q) || out == NULL) return ML_ERR_INVALID_ARG;
    if (!ml_interval_is_supported ()) return ML_ERR_UNSUPPORTED;
    st = mpq_compare_double (q, -DBL_MAX, &cmp);
    if (st != ML_SUCCESS) return st;
    if (cmp < 0) return ML_ERR_OVERFLOW;
    st = mpq_compare_double (q, DBL_MAX, &cmp);
    if (st != ML_SUCCESS) return st;
    if (cmp > 0) return ML_ERR_OVERFLOW;
    lo = double_order_key (-DBL_MAX);
    hi = double_order_key (DBL_MAX);
    while (lo < hi) {
        uint64_t distance = hi - lo;
        uint64_t mid = lo + distance / 2u + distance % 2u;
        d = double_from_order_key (mid);
        st = mpq_compare_double (q, d, &cmp);
        if (st != ML_SUCCESS) return st;
        if (cmp >= 0) lo = mid;
        else hi = mid - 1u;
    }
    d = double_from_order_key (lo);
    st = mpq_compare_double (q, d, &cmp);
    if (st != ML_SUCCESS) return st;
    if (cmp == 0) {
        r.lo = d;
        r.hi = d;
    } else {
        r.lo = d;
        r.hi = nextafter (d, ml_make_inf (0));
        if (!ml_isfinite (r.hi)) return ML_ERR_OVERFLOW;
    }
    *out = r;
    return ML_SUCCESS;
}
