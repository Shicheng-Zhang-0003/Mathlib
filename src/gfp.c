#include "ml_gfp.h"

#include <limits.h>

#define GFP_COOKIE UINT32_C (0x6d4f2a91)
void ml_gfp_ctx_init (ml_gfp_ctx_t *ctx) {
    if (ctx == NULL) return;
    ctx->p = 0u;
    ctx->initialized = 0u;
}
static int gfp_is_prime (uint32_t p) {
    uint32_t d;
    if (p < 2u) return 0;
    if (p == 2u) return 1;
    if ((p & 1u) == 0u) return 0;
    for (d = 3u; (uint64_t) d * d <= p; d += 2u)
        if (p % d == 0u) return 0;
    return 1;
}
ml_status_t ml_gfp_ctx_set_prime (ml_gfp_ctx_t *ctx, uint32_t p) {
    ml_gfp_ctx_t t;
    if (ctx == NULL || !gfp_is_prime (p)) return ML_ERR_INVALID_ARG;
    t.p = p;
    t.initialized = p ^ GFP_COOKIE;
    *ctx = t;
    return ML_SUCCESS;
}
int ml_gfp_ctx_is_valid (const ml_gfp_ctx_t *ctx) {
    return ctx != NULL && ctx->p >= 2u && ctx->initialized == (ctx->p ^ GFP_COOKIE);
}
static int gfp_inputs_valid (const ml_gfp_ctx_t *ctx, uint32_t a,
                             uint32_t *out) {
    return ml_gfp_ctx_is_valid (ctx) && a < ctx->p && out != NULL;
}
static uint32_t gfp_mul_raw (uint32_t p, uint32_t a, uint32_t b) {
    return (uint32_t) (((uint64_t) a * b) % p);
}
static uint32_t gfp_pow_raw (uint32_t p, uint32_t a, uint64_t exponent) {
    uint32_t result = 1u % p;
    while (exponent != 0u) {
        if (exponent & 1u) result = gfp_mul_raw (p, result, a);
        exponent >>= 1;
        if (exponent != 0u) a = gfp_mul_raw (p, a, a);
    }
    return result;
}
ml_status_t ml_gfp_from_i64 (const ml_gfp_ctx_t *ctx, int64_t value,
                             uint32_t *out) {
    uint64_t magnitude, rem;
    uint32_t r;
    if (!ml_gfp_ctx_is_valid (ctx) || out == NULL) return ML_ERR_INVALID_ARG;
    magnitude = value < 0 ? (uint64_t) (-(value + 1)) + 1u : (uint64_t) value;
    rem = magnitude % ctx->p;
    r = (uint32_t) rem;
    if (value < 0 && r != 0u) r = ctx->p - r;
    *out = r;
    return ML_SUCCESS;
}
ml_status_t ml_gfp_add (const ml_gfp_ctx_t *ctx, uint32_t a, uint32_t b,
                        uint32_t *out) {
    uint32_t r;
    if (!gfp_inputs_valid (ctx, a, out) || b >= ctx->p) return ML_ERR_INVALID_ARG;
    r = (uint32_t) (((uint64_t) a + b) % ctx->p);
    *out = r;
    return ML_SUCCESS;
}
ml_status_t ml_gfp_sub (const ml_gfp_ctx_t *ctx, uint32_t a, uint32_t b,
                        uint32_t *out) {
    uint32_t r;
    if (!gfp_inputs_valid (ctx, a, out) || b >= ctx->p) return ML_ERR_INVALID_ARG;
    r = a >= b ? a - b : ctx->p - (b - a);
    *out = r;
    return ML_SUCCESS;
}
ml_status_t ml_gfp_mul (const ml_gfp_ctx_t *ctx, uint32_t a, uint32_t b,
                        uint32_t *out) {
    uint32_t r;
    if (!gfp_inputs_valid (ctx, a, out) || b >= ctx->p) return ML_ERR_INVALID_ARG;
    r = gfp_mul_raw (ctx->p, a, b);
    *out = r;
    return ML_SUCCESS;
}
ml_status_t ml_gfp_pow (const ml_gfp_ctx_t *ctx, uint32_t a,
                        uint64_t exponent, uint32_t *out) {
    uint32_t r;
    if (!gfp_inputs_valid (ctx, a, out)) return ML_ERR_INVALID_ARG;
    r = gfp_pow_raw (ctx->p, a, exponent);
    *out = r;
    return ML_SUCCESS;
}
ml_status_t ml_gfp_inv (const ml_gfp_ctx_t *ctx, uint32_t a, uint32_t *out) {
    uint32_t r;
    if (!gfp_inputs_valid (ctx, a, out)) return ML_ERR_INVALID_ARG;
    if (a == 0u) return ML_ERR_SINGULAR;
    r = gfp_pow_raw (ctx->p, a, (uint64_t) ctx->p - 2u);
    *out = r;
    return ML_SUCCESS;
}
ml_status_t ml_gfp_div (const ml_gfp_ctx_t *ctx, uint32_t a, uint32_t b,
                        uint32_t *out) {
    uint32_t inv, r;
    if (!gfp_inputs_valid (ctx, a, out) || b >= ctx->p) return ML_ERR_INVALID_ARG;
    if (b == 0u) return ML_ERR_SINGULAR;
    inv = gfp_pow_raw (ctx->p, b, (uint64_t) ctx->p - 2u);
    r = gfp_mul_raw (ctx->p, a, inv);
    *out = r;
    return ML_SUCCESS;
}
static ml_status_t gfp_legendre_raw (const ml_gfp_ctx_t *ctx, uint32_t a,
                                     int *symbol) {
    uint32_t r;
    if (!ml_gfp_ctx_is_valid (ctx) || a >= ctx->p || symbol == NULL)
        return ML_ERR_INVALID_ARG;
    if (a == 0u) {
        *symbol = 0;
        return ML_SUCCESS;
    }
    r = gfp_pow_raw (ctx->p, a, ((uint64_t) ctx->p - 1u) / 2u);
    if (r == 1u) *symbol = 1;
    else if (r == ctx->p - 1u) *symbol = -1;
    else return ML_ERR_INTERNAL;
    return ML_SUCCESS;
}
ml_status_t ml_gfp_legendre (const ml_gfp_ctx_t *ctx, uint32_t a, int *symbol) {
    return gfp_legendre_raw (ctx, a, symbol);
}
ml_status_t ml_gfp_sqrt (const ml_gfp_ctx_t *ctx, uint32_t a, uint32_t *root) {
    uint32_t p, x, t, c, z, q;
    unsigned s = 0u, m;
    int symbol;
    ml_status_t st;
    if (!gfp_inputs_valid (ctx, a, root)) return ML_ERR_INVALID_ARG;
    p = ctx->p;
    if (a == 0u || p == 2u) {
        *root = a;
        return ML_SUCCESS;
    }
    st = gfp_legendre_raw (ctx, a, &symbol);
    if (st != ML_SUCCESS) return st;
    if (symbol < 0) return ML_ERR_SINGULAR;
    if ((p & 3u) == 3u) {
        x = gfp_pow_raw (p, a, ((uint64_t) p + 1u) / 4u);
    } else {
        q = p - 1u;
        while ((q & 1u) == 0u) {
            q >>= 1;
            s++;
        }
        z = 2u;
        while (z < p && gfp_pow_raw (p, z, ((uint64_t) p - 1u) / 2u) != p - 1u)
            z++;
        if (z == p) return ML_ERR_INTERNAL;
        x = gfp_pow_raw (p, a, ((uint64_t) q + 1u) / 2u);
        t = gfp_pow_raw (p, a, q);
        c = gfp_pow_raw (p, z, q);
        m = s;
        while (t != 1u) {
            unsigned i;
            uint32_t b, b2;
            uint32_t t2 = t;
            for (i = 1u; i < m; i++) {
                t2 = gfp_mul_raw (p, t2, t2);
                if (t2 == 1u) break;
            }
            if (i == m || m <= i) return ML_ERR_INTERNAL;
            b = gfp_pow_raw (p, c, UINT64_C (1) << (m - i - 1u));
            x = gfp_mul_raw (p, x, b);
            b2 = gfp_mul_raw (p, b, b);
            t = gfp_mul_raw (p, t, b2);
            c = b2;
            m = i;
        }
    }
    if (x > p - x) x = p - x;
    *root = x;
    return ML_SUCCESS;
}
