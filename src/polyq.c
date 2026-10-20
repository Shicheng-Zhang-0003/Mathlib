#include "ml_polyq.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
void ml_polyq_init (ml_polyq_t *p) {
    if (p == NULL) return;
    p->coefficient = NULL;
    p->size = 0u;
    p->capacity = 0u;
}
void ml_polyq_clear (ml_polyq_t *p) {
    size_t i;
    if (p == NULL) return;
    for (i = 0u; i < p->capacity; i++) ml_mpq_clear (&p->coefficient [i]);
    free (p->coefficient);
    ml_polyq_init (p);
}
int ml_polyq_is_valid (const ml_polyq_t *p) {
    size_t i;
    if (p == NULL || p->size > p->capacity ||
        ((p->capacity == 0u) != (p->coefficient == NULL))) return 0;
    for (i = 0u; i < p->size; i++) {
        if (!ml_mpq_is_valid (&p->coefficient [i])) return 0;
    }
    if (p->size != 0u && p->coefficient [p->size - 1u].numerator.size == 0u)
        return 0;
    return 1;
}
static int polyq_initialized (const ml_polyq_t *p) {
    return p != NULL && p->size == 0u && p->capacity == 0u &&
           p->coefficient == NULL;
}
int ml_polyq_is_zero (const ml_polyq_t *p) {
    return ml_polyq_is_valid (p) && p->size == 0u;
}
ptrdiff_t ml_polyq_degree (const ml_polyq_t *p) {
    if (!ml_polyq_is_valid (p)) return (ptrdiff_t) -2;
    return p->size == 0u ? (ptrdiff_t) -1 : (ptrdiff_t) (p->size - 1u);
}
static void polyq_commit (ml_polyq_t *out, ml_polyq_t *tmp) {
    ml_polyq_clear (out);
    *out = *tmp;
    ml_polyq_init (tmp);
}
static ml_status_t polyq_reserve (ml_polyq_t *p, size_t need) {
    ml_mpq_t *items;
    size_t cap, i, initialized = 0u;
    ml_status_t st;
    if (need <= p->capacity) return ML_SUCCESS;
    if (need > SIZE_MAX / sizeof (ml_mpq_t)) return ML_ERR_OVERFLOW;
    cap = p->capacity ? p->capacity : 4u;
    while (cap < need) {
        if (cap > SIZE_MAX / 2u) {
            cap = need;
            break;
        }
        cap *= 2u;
    }
    if (cap > SIZE_MAX / sizeof (ml_mpq_t)) cap = need;
    items = (ml_mpq_t *) calloc (cap, sizeof (ml_mpq_t));
    if (items == NULL) return ML_ERR_WORKSPACE;
    if (p->capacity != 0u)
        memcpy (items, p->coefficient, p->capacity * sizeof (ml_mpq_t));
    for (i = p->capacity; i < cap; i++) {
        ml_mpq_init (&items [i]);
        initialized++;
        st = ml_mpq_set_i64 (&items [i], 0);
        if (st != ML_SUCCESS) {
            size_t j;
            for (j = p->capacity; j < p->capacity + initialized; j++)
                ml_mpq_clear (&items [j]);
            free (items);
            return st;
        }
    }
    free (p->coefficient); /* ownership of active coefficient values moved */
    p->coefficient = items;
    p->capacity = cap;
    return ML_SUCCESS;
}
static void polyq_trim (ml_polyq_t *p) {
    while (p->size != 0u && p->coefficient [p->size - 1u].numerator.size == 0u)
        p->size--;
}
ml_status_t ml_polyq_copy (ml_polyq_t *r, const ml_polyq_t *p) {
    ml_polyq_t t;
    size_t i;
    ml_status_t st;
    if (!polyq_initialized (r) && !ml_polyq_is_valid (r)) return ML_ERR_INVALID_ARG;
    if (!ml_polyq_is_valid (p)) return ML_ERR_INVALID_ARG;
    if (r == p) return ML_SUCCESS;
    ml_polyq_init (&t);
    st = polyq_reserve (&t, p->size);
    if (st != ML_SUCCESS) goto done;
    for (i = 0u; i < p->size; i++) {
        st = ml_mpq_copy (&t.coefficient [i], &p->coefficient [i]);
        if (st != ML_SUCCESS) goto done;
    }
    t.size = p->size;
    polyq_commit (r, &t);
done:
    ml_polyq_clear (&t);
    return st;
}
ml_status_t ml_polyq_set_coefficient (ml_polyq_t *p, size_t degree,
                                      const ml_mpq_t *value) {
    ml_mpq_t t;
    ml_status_t st;
    if (!ml_polyq_is_valid (p) || !ml_mpq_is_valid (value)) return ML_ERR_INVALID_ARG;
    if (degree == SIZE_MAX) return ML_ERR_OVERFLOW;
    ml_mpq_init (&t);
    st = ml_mpq_set_i64 (&t, 0);
    if (st == ML_SUCCESS) st = ml_mpq_copy (&t, value);
    if (st != ML_SUCCESS) {
        ml_mpq_clear (&t);
        return st;
    }
    st = polyq_reserve (p, degree + 1u);
    if (st == ML_SUCCESS) st = ml_mpq_copy (&p->coefficient [degree], &t);
    if (st == ML_SUCCESS) {
        if (p->size <= degree) p->size = degree + 1u;
        polyq_trim (p);
    }
    ml_mpq_clear (&t);
    return st;
}
ml_status_t ml_polyq_get_coefficient (const ml_polyq_t *p, size_t degree,
                                      ml_mpq_t *value) {
    if (!ml_polyq_is_valid (p) || value == NULL ||
        (!ml_mpq_is_valid (value) &&
         (!ml_mpz_is_valid (&value->numerator) || !ml_mpz_is_valid (&value->denominator))))
        return ML_ERR_INVALID_ARG;
    if (degree >= p->size) return ml_mpq_set_i64 (value, 0);
    return ml_mpq_copy (value, &p->coefficient [degree]);
}
ml_status_t ml_polyq_equal (const ml_polyq_t *a, const ml_polyq_t *b,
                            int *equal) {
    size_t i;
    int cmp;
    ml_status_t st;
    if (!ml_polyq_is_valid (a) || !ml_polyq_is_valid (b) || equal == NULL)
        return ML_ERR_INVALID_ARG;
    if (a->size != b->size) {
        *equal = 0;
        return ML_SUCCESS;
    }
    for (i = 0u; i < a->size; i++) {
        st = ml_mpq_compare (&a->coefficient [i], &b->coefficient [i], &cmp);
        if (st != ML_SUCCESS) return st;
        if (cmp != 0) {
            *equal = 0;
            return ML_SUCCESS;
        }
    }
    *equal = 1;
    return ML_SUCCESS;
}
static ml_status_t polyq_binary (ml_polyq_t *r, const ml_polyq_t *a,
                                 const ml_polyq_t *b, int subtract) {
    ml_polyq_t t;
    ml_mpq_t zero;
    size_t n, i;
    ml_status_t st;
    if ((!polyq_initialized (r) && !ml_polyq_is_valid (r)) ||
        !ml_polyq_is_valid (a) || !ml_polyq_is_valid (b)) return ML_ERR_INVALID_ARG;
    ml_polyq_init (&t);
    ml_mpq_init (&zero);
    st = ml_mpq_set_i64 (&zero, 0);
    if (st != ML_SUCCESS) goto done;
    n = a->size > b->size ? a->size : b->size;
    st = polyq_reserve (&t, n);
    if (st != ML_SUCCESS) goto done;
    for (i = 0u; i < n; i++) {
        const ml_mpq_t *av = i < a->size ? &a->coefficient [i] : &zero;
        const ml_mpq_t *bv = i < b->size ? &b->coefficient [i] : &zero;
        st = subtract ? ml_mpq_sub (&t.coefficient [i], av, bv)
                      : ml_mpq_add (&t.coefficient [i], av, bv);
        if (st != ML_SUCCESS) goto done;
    }
    t.size = n;
    polyq_trim (&t);
    polyq_commit (r, &t);
done:
    ml_polyq_clear (&t);
    ml_mpq_clear (&zero);
    return st;
}
ml_status_t ml_polyq_add (ml_polyq_t *r, const ml_polyq_t *a,
                          const ml_polyq_t *b) {
    return polyq_binary (r, a, b, 0);
}
ml_status_t ml_polyq_sub (ml_polyq_t *r, const ml_polyq_t *a,
                          const ml_polyq_t *b) {
    return polyq_binary (r, a, b, 1);
}
ml_status_t ml_polyq_mul (ml_polyq_t *r, const ml_polyq_t *a,
                          const ml_polyq_t *b) {
    ml_polyq_t t;
    ml_mpq_t product;
    size_t n, i, j;
    ml_status_t st = ML_SUCCESS;
    if ((!polyq_initialized (r) && !ml_polyq_is_valid (r)) ||
        !ml_polyq_is_valid (a) || !ml_polyq_is_valid (b)) return ML_ERR_INVALID_ARG;
    ml_polyq_init (&t);
    ml_mpq_init (&product);
    st = ml_mpq_set_i64 (&product, 0);
    if (st != ML_SUCCESS) goto done;
    if (a->size == 0u || b->size == 0u) {
        polyq_commit (r, &t);
        goto done;
    }
    if (a->size > SIZE_MAX - b->size + 1u) {
        st = ML_ERR_OVERFLOW;
        goto done;
    }
    n = a->size + b->size - 1u;
    st = polyq_reserve (&t, n);
    if (st != ML_SUCCESS) goto done;
    t.size = n;
    for (i = 0u; i < a->size; i++) {
        for (j = 0u; j < b->size; j++) {
            size_t k = i + j;
            st = ml_mpq_mul (&product, &a->coefficient [i], &b->coefficient [j]);
            if (st != ML_SUCCESS) goto done;
            st = ml_mpq_add (&t.coefficient [k], &t.coefficient [k], &product);
            if (st != ML_SUCCESS) goto done;
        }
    }
    polyq_trim (&t);
    polyq_commit (r, &t);
done:
    ml_polyq_clear (&t);
    ml_mpq_clear (&product);
    return st;
}
ml_status_t ml_polyq_derivative (ml_polyq_t *r, const ml_polyq_t *p) {
    ml_polyq_t t;
    ml_mpq_t factor;
    size_t i, n;
    ml_status_t st = ML_SUCCESS;
    if ((!polyq_initialized (r) && !ml_polyq_is_valid (r)) || !ml_polyq_is_valid (p))
        return ML_ERR_INVALID_ARG;
    ml_polyq_init (&t);
    ml_mpq_init (&factor);
    st = ml_mpq_set_i64 (&factor, 0);
    if (st != ML_SUCCESS) goto done;
    n = p->size > 0u ? p->size - 1u : 0u;
    st = polyq_reserve (&t, n);
    if (st != ML_SUCCESS) goto done;
    for (i = 1u; i < p->size; i++) {
        if (i > (size_t) INT64_MAX) {
            st = ML_ERR_OVERFLOW;
            goto done;
        }
        st = ml_mpq_set_i64 (&factor, (int64_t) i);
        if (st != ML_SUCCESS) goto done;
        st = ml_mpq_mul (&t.coefficient [i - 1u], &p->coefficient [i], &factor);
        if (st != ML_SUCCESS) goto done;
    }
    t.size = n;
    polyq_trim (&t);
    polyq_commit (r, &t);
done:
    ml_polyq_clear (&t);
    ml_mpq_clear (&factor);
    return st;
}
ml_status_t ml_polyq_evaluate (ml_mpq_t *value, const ml_polyq_t *p,
                               const ml_mpq_t *x) {
    ml_mpq_t acc;
    size_t i;
    ml_status_t st;
    if (value == NULL || !ml_mpz_is_valid (&value->numerator) ||
        !ml_mpz_is_valid (&value->denominator) || !ml_polyq_is_valid (p) ||
        !ml_mpq_is_valid (x)) return ML_ERR_INVALID_ARG;
    ml_mpq_init (&acc);
    st = ml_mpq_set_i64 (&acc, 0);
    i = p->size;
    while (st == ML_SUCCESS && i > 0u) {
        i--;
        st = ml_mpq_mul (&acc, &acc, x);
        if (st == ML_SUCCESS)
            st = ml_mpq_add (&acc, &acc, &p->coefficient [i]);
    }
    if (st == ML_SUCCESS) st = ml_mpq_copy (value, &acc);
    ml_mpq_clear (&acc);
    return st;
}
ml_status_t ml_polyq_divrem (ml_polyq_t *quotient, ml_polyq_t *remainder,
                             const ml_polyq_t *dividend,
                             const ml_polyq_t *divisor) {
    ml_polyq_t q, rem;
    ml_mpq_t factor, product;
    size_t i, old_size;
    ml_status_t st = ML_SUCCESS;
    if (quotient == remainder ||
        ((!polyq_initialized (quotient) && !ml_polyq_is_valid (quotient))) ||
        ((!polyq_initialized (remainder) && !ml_polyq_is_valid (remainder))) ||
        !ml_polyq_is_valid (dividend) || !ml_polyq_is_valid (divisor))
        return ML_ERR_INVALID_ARG;
    if (divisor->size == 0u) return ML_ERR_SINGULAR;
    ml_polyq_init (&q);
    ml_polyq_init (&rem);
    ml_mpq_init (&factor);
    ml_mpq_init (&product);
    st = ml_mpq_set_i64 (&factor, 0);
    if (st == ML_SUCCESS) st = ml_mpq_set_i64 (&product, 0);
    if (st == ML_SUCCESS) st = ml_polyq_copy (&rem, dividend);
    if (st != ML_SUCCESS) goto done;
    if (rem.size >= divisor->size) {
        size_t qn = rem.size - divisor->size + 1u;
        st = polyq_reserve (&q, qn);
        if (st != ML_SUCCESS) goto done;
        q.size = qn;
    }
    while (rem.size >= divisor->size && rem.size != 0u) {
        size_t k = rem.size - divisor->size;
        old_size = rem.size;
        st = ml_mpq_div (&factor, &rem.coefficient [rem.size - 1u],
                         &divisor->coefficient [divisor->size - 1u]);
        if (st != ML_SUCCESS) goto done;
        st = ml_mpq_copy (&q.coefficient [k], &factor);
        if (st != ML_SUCCESS) goto done;
        for (i = 0u; i < divisor->size; i++) {
            st = ml_mpq_mul (&product, &factor, &divisor->coefficient [i]);
            if (st != ML_SUCCESS) goto done;
            st = ml_mpq_sub (&rem.coefficient [k + i], &rem.coefficient [k + i],
                             &product);
            if (st != ML_SUCCESS) goto done;
        }
        polyq_trim (&rem);
        if (rem.size >= old_size) {
            st = ML_ERR_INTERNAL;
            goto done;
        }
    }
    polyq_trim (&q);
    polyq_commit (quotient, &q);
    polyq_commit (remainder, &rem);
done:
    ml_polyq_clear (&q);
    ml_polyq_clear (&rem);
    ml_mpq_clear (&factor);
    ml_mpq_clear (&product);
    return st;
}
ml_status_t ml_polyq_gcd (ml_polyq_t *gcd, const ml_polyq_t *a,
                          const ml_polyq_t *b) {
    ml_polyq_t x, y, q, rem;
    ml_mpq_t lead;
    ml_status_t st;
    size_t i;
    if ((!polyq_initialized (gcd) && !ml_polyq_is_valid (gcd)) ||
        !ml_polyq_is_valid (a) || !ml_polyq_is_valid (b)) return ML_ERR_INVALID_ARG;
    ml_polyq_init (&x);
    ml_polyq_init (&y);
    ml_polyq_init (&q);
    ml_polyq_init (&rem);
    ml_mpq_init (&lead);
    st = ml_polyq_copy (&x, a);
    if (st == ML_SUCCESS) st = ml_polyq_copy (&y, b);
    while (st == ML_SUCCESS && y.size != 0u) {
        st = ml_polyq_divrem (&q, &rem, &x, &y);
        if (st != ML_SUCCESS) break;
        {
            ml_polyq_t tmp = x;
            x = y;
            y = rem;
            rem = tmp;
        }
        ml_polyq_clear (&q);
    }
    if (st != ML_SUCCESS) goto done;
    if (x.size != 0u) {
        st = ml_mpq_copy (&lead, &x.coefficient [x.size - 1u]);
        if (st != ML_SUCCESS) goto done;
        for (i = 0u; i < x.size; i++) {
            st = ml_mpq_div (&x.coefficient [i], &x.coefficient [i], &lead);
            if (st != ML_SUCCESS) goto done;
        }
    }
    polyq_commit (gcd, &x);
done:
    ml_polyq_clear (&x);
    ml_polyq_clear (&y);
    ml_polyq_clear (&q);
    ml_polyq_clear (&rem);
    ml_mpq_clear (&lead);
    return st;
}
ml_status_t ml_polyq_verify_division (const ml_polyq_t *dividend,
                                      const ml_polyq_t *divisor,
                                      const ml_polyq_t *quotient,
                                      const ml_polyq_t *remainder,
                                      int *valid) {
    ml_polyq_t product, sum;
    ml_status_t st;
    int equal = 0;
    if (valid == NULL) return ML_ERR_INVALID_ARG;
    *valid = 0;
    if (!ml_polyq_is_valid (dividend) || !ml_polyq_is_valid (divisor) ||
        !ml_polyq_is_valid (quotient) || !ml_polyq_is_valid (remainder))
        return ML_ERR_INVALID_ARG;
    if (divisor->size == 0u) return ML_ERR_SINGULAR;
    if (remainder->size >= divisor->size) return ML_SUCCESS;
    ml_polyq_init (&product);
    ml_polyq_init (&sum);
    st = ml_polyq_mul (&product, divisor, quotient);
    if (st == ML_SUCCESS) st = ml_polyq_add (&sum, &product, remainder);
    if (st == ML_SUCCESS) st = ml_polyq_equal (&sum, dividend, &equal);
    if (st == ML_SUCCESS) *valid = equal;
    ml_polyq_clear (&product);
    ml_polyq_clear (&sum);
    return st;
}
