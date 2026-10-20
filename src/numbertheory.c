#include "ml_compiler.h"
#include "ml_numbertheory.h"
#include "ml_combinatorics.h"
#include <stdlib.h>
static uint64_t ml_abs_u64 (int64_t v) {
    if (v >= 0) return (uint64_t) v;
    if (v == INT64_MIN) return (uint64_t) INT64_MAX + 1ULL;
    return (uint64_t) (-v);
}
ML_API int ml_legendre_symbol (int64_t a, uint64_t p) {
    if (p < 3 || (p & 1ULL) == 0ULL) return 0;
    if (!ml_is_prime64 (p)) return 0;
    {
        int64_t pp = (int64_t) p;
        int64_t t = a % pp;
        if (t < 0) t += pp;
        uint64_t am = (uint64_t) t;
        if (am == 0) return 0;
        uint64_t r = ml_modpow (am, (p - 1) / 2, p);
        if (r == 1) return 1;
        if (r == p - 1) return -1;
        return 0;
    }
}
ML_API int ml_jacobi_symbol (int64_t a, uint64_t n) {
    if (n == 0 || (n & 1ULL) == 0ULL) return 0;
    if (n == 1) return 1;
    {
        uint64_t aa, nn = n;
        int t = 1;
        /* (a/n) depends only on a mod n, reduced into [0,n).  Reducing in
         * the magnitude domain keeps this exact for a = INT64_MIN and for
         * n > 2^63-1, and the loop below stays entirely unsigned. */
        if (a < 0) {
            uint64_t m = (uint64_t) (-(a + 1)) + 1ULL; /* |a|, no INT64_MIN UB */
            uint64_t r = m % nn;
            aa = (r == 0ULL) ? 0ULL : nn - r;
        } else {
            aa = (uint64_t) a % nn;
        }
        while (aa != 0) {
            while ((aa & 1ULL) == 0ULL) {
                uint64_t r;
                aa >>= 1;
                r = nn & 7ULL;
                if (r == 3ULL || r == 5ULL) t = -t;
            }
            {
                uint64_t tmp = aa;
                aa = nn;
                nn = tmp;
                if ((aa & 3ULL) == 3ULL && (nn & 3ULL) == 3ULL) t = -t;
                aa %= nn;
            }
        }
        return (nn == 1) ? t : 0;
    }
}
ML_API int ml_kronecker_symbol (int64_t a, int64_t n) {
    if (n == 0) {
        uint64_t aa = ml_abs_u64 (a);
        return (aa == 1) ? 1 : 0;
    }
    if (n == 1) return 1;
    if (n == -1) return (a < 0) ? -1 : 1;
    {
        int t = 1;
        /* DESPOT-AUDIT: INT64_MIN cannot be negated (UB). Map magnitude via
         * ml_abs_u64 first; sign handling for n<0 uses the same rule. */
        uint64_t nn;
        if (n < 0) {
            nn = ml_abs_u64 (n);
            if (a < 0) t = -t;
        } else {
            nn = (uint64_t) n;
        }
        /* Factor out powers of 2: (a/2) = 0 if even, else +-1 by a mod 8. */
        int e = 0;
        while ((nn & 1ULL) == 0ULL) {
            nn >>= 1;
            e++;
        }
        if (e > 0) {
            uint64_t am8 = ml_abs_u64 (a) & 7ULL;
            if ((a & 1LL) == 0) return 0;
            int a2 = 0;
            if (am8 == 1 || am8 == 7) a2 = 1;
            else if (am8 == 3 || am8 == 5) a2 = -1;
            if (e % 2 == 1) t *= a2;
            else if (a2 == 0) return 0;
        }
        if (nn == 1) return t;
        return t * ml_jacobi_symbol (a, nn);
    }
}
ML_API int ml_mobius (int64_t n) {
    /* Contract: returns -1/0/1. Returns 2 as UNRESOLVED sentinel when the
     * cofactor beyond the 1e6 trial bound is composite with unknown
     * square-freeness (documented in ml_numbertheory.h). Callers must treat
     * 2 as "unknown", never as squarefree/non-squarefree. */
    if (n == 0) return 0;
    uint64_t m = ml_abs_u64 (n);
    if (m == 1) return 1;
    {
        int prime_factors = 0;
        uint64_t tmp = m;
        for (uint64_t p = 2; p <= tmp / p; p += (p == 2 ? 1 : 2)) {
            if (p > 1000000ULL) {
                /* Too large for trial: use prime test to decide. If the
                 * cofactor is prime, count it; else unknown -> 2. */
                if (ml_is_prime64 (tmp)) {
                    prime_factors++;
                    tmp = 1;
                    break;
                }
                return 2;
            }
            if (tmp % p == 0) {
                tmp /= p;
                if (tmp % p == 0) return 0;
                prime_factors++;
            }
        }
        if (tmp > 1) prime_factors++;
        return (prime_factors % 2 == 0) ? 1 : -1;
    }
}
ML_API uint64_t ml_mult_order (uint64_t a, uint64_t m) {
    if (m <= 1) return 0;
    a %= m;
    if (ml_gcd (a, m) != 1) return 0;
    {
        uint64_t phi = ml_phi (m);
        if (phi == 0) return 0;
        /* Order divides phi: strip prime factors. Trial to 1e6 covers the
         * small-factor path; any remaining cofactor tmp>1 (prime or
         * composite) is stripped one prime at a time via gcd-division:
         * while ord % tmp == 0 and a^(ord/tmp) == 1, divide. If tmp is
         * composite this still converges to the true order because each
         * successful division is re-tested until no prime factor of tmp
         * divides the cofinal ord. Residual composite tmp that is NOT a
         * divisor of ord is ignored (cannot divide the order). */
        uint64_t ord = phi, tmp = phi;
        for (uint64_t p = 2; p <= tmp / p; p += (p == 2 ? 1 : 2)) {
            if (p > 1000000ULL) break;
            if (tmp % p == 0) {
                while (tmp % p == 0) tmp /= p;
                while (ord % p == 0 && ml_modpow (a, ord / p, m) == 1) ord /= p;
            }
        }
        while (tmp > 1 && ord % tmp == 0 && ml_modpow (a, ord / tmp, m) == 1) {
            /* tmp may be composite: factor one step via gcd with ord/tmp
             * probe. If a^(ord/tmp) == 1, ord/tmp is still a multiple of
             * the order, so accept and re-derive tmp = gcd(tmp, ord). */
            ord /= tmp;
            /* Recompute residual: tmp := gcd(tmp, ord) loop guard so a
             * composite tmp with a prime factor already removed cannot
             * over-divide. */
            uint64_t g = ml_gcd (tmp, ord);
            if (g <= 1) break;
            tmp = g;
            if (tmp == 1) break;
            /* Continue stripping the reduced tmp by trial below. */
            for (uint64_t p = 2; p <= tmp / p && p <= 1000000ULL; p += (p == 2 ? 1 : 2)) {
                if (tmp % p == 0) {
                    while (tmp % p == 0) tmp /= p;
                    while (ord % p == 0 && ml_modpow (a, ord / p, m) == 1) ord /= p;
                }
            }
        }
        return ord;
    }
}
ML_API int ml_is_primitive_root (uint64_t a, uint64_t p) {
    if (p < 2 || !ml_is_prime64 (p)) return 0;
    a %= p;
    if (a == 0) return 0;
    return ml_mult_order (a, p) == p - 1 ? 1 : 0;
}
ML_API ml_status_t ml_crt2 (uint64_t a1, uint64_t m1, uint64_t a2, uint64_t m2,
                            uint64_t *x, uint64_t *lcm) {
    if (ML_UNLIKELY (!x || !lcm || m1 == 0 || m2 == 0)) return ML_ERR_INVALID_ARG;
    {
        uint64_t g = ml_gcd (m1, m2);
        uint64_t d1 = (a2 >= a1) ? a2 - a1 : a1 - a2;
        if (d1 % g != 0) return ML_ERR_SINGULAR;
        {
            uint64_t l = m1 / g * m2;
            if (l / m2 != m1 / g) return ML_ERR_WORKSPACE;
            /* Solve m1*t = a2-a1 (mod m2): reduce by g, invert.
             * DESPOT-AUDIT: diff must be reduced mod q=m2/g, not mod m2.
             * Previous (m2-((a1-a2)/g)%m2)%m2 was wrong for g>1,a2<a1
             * (e.g. x=0 mod6, x=4 mod10 gave 6 instead of 24 mod30). */
            uint64_t p = m1 / g, q = m2 / g;
            uint64_t diff = (a2 >= a1) ? (a2 - a1) / g : (q - (((a1 - a2) / g) % q)) % q;
            {
                uint64_t inv = ml_modinv (p % q, q);
                if (inv == UINT64_MAX) {
                    /* Fallback extended Euclid on reduced moduli.
                     * egcd is int64-only: reduced moduli beyond INT64_MAX
                     * cannot use this path (would wrap). modinv already
                     * handles the full u64 range, so reaching here with
                     * huge q means genuinely non-invertible. */
                    if (p > (uint64_t) INT64_MAX || q > (uint64_t) INT64_MAX)
                        return ML_ERR_SINGULAR;
                    int64_t ex = 0, ey = 0;
                    int64_t gg = ml_egcd ((int64_t) p, (int64_t) q, &ex, &ey);
                    if (gg != 1 && gg != -1) return ML_ERR_SINGULAR;
                    ex %= (int64_t) q;
                    if (ex < 0) ex += (int64_t) q;
                    inv = (uint64_t) ex;
                }
                __uint128_t t = (__uint128_t) diff * inv % q;
                __uint128_t r = (__uint128_t) a1 + (__uint128_t) m1 * t;
                *x = (uint64_t) (r % l);
                *lcm = l;
                return ML_SUCCESS;
            }
        }
    }
}
ML_API uint64_t ml_prime_pi (uint64_t n) {
    if (n < 2) return 0;
    if (n > 100000000ULL) return UINT64_MAX;
    {
        uint64_t nn = n;
        /* Heap-backed odd-only sieve sized by the limit: avoids a 50MB
         * static BSS image and is thread-safe (no shared mutable state).
         * DESPOT-FIX: size=(nn+1)/2 so odd nn itself is sieved/counted.
         * Old nn/2 missed nn when nn was an odd prime (pi(3)=1, pi(5)=2). */
        uint64_t size = (nn + 1) / 2;
        /* R5: SIZE_MAX guard (50MB max by the 1e8 cap, so defense-in-depth
         * on 32-bit). OOM maps to UINT64_MAX sentinel (documented: caller
         * must treat UINT64_MAX as "unknown", not a count). */
        if ((size_t) size != size || size > (size_t) -1 / sizeof (uint8_t)) return UINT64_MAX;
        uint8_t *bits = (uint8_t *) malloc ((size_t) size * sizeof (uint8_t));
        if (!bits) return UINT64_MAX;
        for (uint64_t i = 0; i < size; i++) bits [i] = 1;
        for (uint64_t p = 3; p * p <= nn; p += 2) {
            if (bits [p / 2]) {
                for (uint64_t q2 = p * p; q2 <= nn; q2 += 2 * p) bits [q2 / 2] = 0;
            }
        }
        uint64_t c = 1;
        for (uint64_t i = 1; i < size; i++) {
            if (bits [i]) c++;
        }
        free (bits);
        return c;
    }
}
ML_API int ml_quadratic_reciprocity_check (uint64_t p, uint64_t q) {
    if (!ml_is_prime64 (p) || !ml_is_prime64 (q) || p < 3 || q < 3 || p == q) return 0;
    {
        int64_t lp = ml_legendre_symbol ((int64_t) p, q);
        int64_t lq = ml_legendre_symbol ((int64_t) q, p);
        int sign = (((p - 1) / 2) & ((q - 1) / 2) & 1) ? -1 : 1;
        return (lp == sign * lq) ? 1 : 0;
    }
}
/* ---- Olympiad kit: LTE, floor-sum, Legendre-n!, continued fractions, Pell.
 * Exact u64 paths with UINT64_MAX sentinel (house convention). */

static uint64_t mln_vp_u64 (uint64_t v, uint64_t p) {
    uint64_t c = 0;
    if (v == 0 || p < 2) return 0;
    while (v % p == 0) {
        v /= p;
        c++;
    }
    return c;
}
static uint64_t mln_isqrt_u64 (uint64_t x) {
    uint64_t r;
    uint64_t y;
    if (x <= 1) return x;
    r = x;
    y = (r + 1) / 2;
    while (y < r) {
        r = y;
        y = (r + x / r) / 2;
    }
    return r;
}
ML_API uint64_t ml_vp_factorial (long long n, uint64_t p) {
    /* Legendre: v_p(n!) = n/p + n/p^2 + ... (prime p only). */
    uint64_t s = 0;
    uint64_t t;
    if (n < 0 || n > 9007199254740992LL) return UINT64_MAX;
    if (p < 2 || !ml_is_prime64 (p)) return UINT64_MAX;
    t = (uint64_t) n;
    while (t != 0) {
        t /= p;
        if (s > UINT64_MAX - t) return UINT64_MAX;
        s += t;
    }
    return s;
}
ML_API uint64_t ml_floor_sum_u64 (uint64_t n, uint64_t m, uint64_t a, uint64_t b) {
    /* AtCoder floor_sum: sum_{0<=i<n} floor((a*i+b)/m). m = 0 invalid.
     * Uses __uint128_t internally (same-file precedent: ml_crt2). */
    __uint128_t ans = 0;
    if (m == 0) return UINT64_MAX;
    while (1) {
        if (a >= m) {
            __uint128_t q = (__uint128_t) (a / m);
            ans += q * (__uint128_t) n * (__uint128_t) (n - 1) / (__uint128_t) 2;
            a %= m;
        }
        if (b >= m) {
            ans += (__uint128_t) n * (__uint128_t) (b / m);
            b %= m;
        }
        {
            __uint128_t y_max = (__uint128_t) a * (__uint128_t) n + (__uint128_t) b;
            if (y_max < (__uint128_t) m) break;
            n = (uint64_t) (y_max / (__uint128_t) m);
            b = (uint64_t) (y_max % (__uint128_t) m);
            {
                uint64_t tmp = m;
                m = a;
                a = tmp;
            }
        }
        if (ans > (__uint128_t) UINT64_MAX) return UINT64_MAX;
    }
    if (ans > (__uint128_t) UINT64_MAX) return UINT64_MAX;
    return (uint64_t) ans;
}
ML_API uint64_t ml_lte_vp (uint64_t a, uint64_t b, uint64_t n, uint64_t p) {
    /* LTE: v_p(a^n-b^n). Odd p: p | a-b, p coprime ab. p = 2: a,b odd;
     * n odd -> v2(a-b); n even -> v2(a-b)+v2(a+b)+v2(n)-1. */
    uint64_t d;
    if (n < 1 || p < 2 || !ml_is_prime64 (p)) return UINT64_MAX;
    if (a == b) return UINT64_MAX;
    d = (a > b) ? a - b : b - a;
    if (d % p != 0) return UINT64_MAX;
    if (a % p == 0) return UINT64_MAX;
    if (p == 2) {
        if ((a & 1ULL) == 0 || (b & 1ULL) == 0) return UINT64_MAX;
        if ((n & 1ULL) != 0) return mln_vp_u64 (d, 2);
        {
            uint64_t s = a + b;
            uint64_t vn;
            if (s < a) return UINT64_MAX;
            vn = mln_vp_u64 (n, 2);
            if (vn == 0) return UINT64_MAX;
            if (mln_vp_u64 (d, 2) > UINT64_MAX - mln_vp_u64 (s, 2)) return UINT64_MAX;
            {
                uint64_t t = mln_vp_u64 (d, 2) + mln_vp_u64 (s, 2);
                if (t > UINT64_MAX - vn) return UINT64_MAX;
                t += vn;
                if (t < 1) return UINT64_MAX;
                return t - 1;
            }
        }
    }
    {
        uint64_t vd = mln_vp_u64 (d, p);
        uint64_t vn = mln_vp_u64 (n, p);
        if (vd > UINT64_MAX - vn) return UINT64_MAX;
        return vd + vn;
    }
}
ML_API int ml_sqrt_cf_period (uint64_t D) {
    /* Period length of sqrt(D) CF. 0 = perfect square, -1 = invalid (D<2),
     * -2 = period exceeds the 10000 bound. */
    uint64_t a0;
    uint64_t m = 0;
    uint64_t d = 1;
    uint64_t a;
    int len = 0;
    if (D < 2) return -1;
    a0 = mln_isqrt_u64 (D);
    if (a0 * a0 == D) return 0;
    a = a0;
    do {
        uint64_t m2;
        uint64_t d2;
        if (a > UINT64_MAX - m) return -2;
        m2 = d * a - m;
        if (m2 * m2 > D) return -2;
        if (D - m2 * m2 == 0) return -2;
        d2 = (D - m2 * m2) / d;
        if (d2 == 0) return -2;
        a = (a0 + m2) / d2;
        m = m2;
        d = d2;
        len++;
        if (len > 10000) return -2;
    } while (a != 2 * a0);
    return len;
}
ML_API ml_status_t ml_cf_convergent_sqrt (uint64_t D, int k, uint64_t *num, uint64_t *den) {
    /* k-th convergent (0-indexed) of sqrt(D) via the CF generator. */
    uint64_t a0;
    uint64_t m = 0;
    uint64_t d = 1;
    uint64_t a;
    uint64_t p0 = 0, p1 = 1;
    uint64_t q0 = 1, q1 = 0;
    int i;
    if (ML_UNLIKELY (!num || !den)) return ML_ERR_INVALID_ARG;
    if (D < 2 || k < 0 || k > 10000) return ML_ERR_INVALID_ARG;
    a0 = mln_isqrt_u64 (D);
    if (a0 * a0 == D) return ML_ERR_INVALID_ARG;
    a = a0;
    for (i = 0; i <= k; i++) {
        uint64_t p2;
        uint64_t q2;
        uint64_t m2;
        uint64_t d2;
        uint64_t anext;
        if (a != 0 && p1 > UINT64_MAX / a) return ML_ERR_OVERFLOW;
        p2 = a * p1 + p0;
        if (p2 < p0) return ML_ERR_OVERFLOW;
        if (a != 0 && q1 > UINT64_MAX / a) return ML_ERR_OVERFLOW;
        q2 = a * q1 + q0;
        if (q2 < q0) return ML_ERR_OVERFLOW;
        p0 = p1;
        p1 = p2;
        q0 = q1;
        q1 = q2;
        if (i == k) break;
        if (a > UINT64_MAX - m) return ML_ERR_OVERFLOW;
        m2 = d * a - m;
        if (m2 * m2 > D) return ML_ERR_OVERFLOW;
        d2 = (D - m2 * m2) / d;
        if (d2 == 0) return ML_ERR_INTERNAL;
        anext = (a0 + m2) / d2;
        m = m2;
        d = d2;
        a = anext;
    }
    *num = p1;
    *den = q1;
    return ML_SUCCESS;
}
ML_API ml_status_t ml_pell_fund (uint64_t D, uint64_t *x, uint64_t *y) {
    /* Minimal x^2-Dy^2 = 1 via convergents (first hit wins; theory bounds
     * it by twice the period). Square D: INVALID_ARG (no solution). */
    uint64_t a0;
    uint64_t m = 0;
    uint64_t d = 1;
    uint64_t a;
    uint64_t p0 = 0, p1 = 1;
    uint64_t q0 = 1, q1 = 0;
    int i;
    if (ML_UNLIKELY (!x || !y)) return ML_ERR_INVALID_ARG;
    if (D < 2) return ML_ERR_INVALID_ARG;
    a0 = mln_isqrt_u64 (D);
    if (a0 * a0 == D) return ML_ERR_INVALID_ARG;
    a = a0;
    for (i = 0; i < 20000; i++) {
        uint64_t p2;
        uint64_t q2;
        uint64_t m2;
        uint64_t d2;
        uint64_t lhs;
        uint64_t rhs;
        if (a != 0 && p1 > UINT64_MAX / a) return ML_ERR_OVERFLOW;
        p2 = a * p1 + p0;
        if (p2 < p0) return ML_ERR_OVERFLOW;
        if (a != 0 && q1 > UINT64_MAX / a) return ML_ERR_OVERFLOW;
        q2 = a * q1 + q0;
        if (q2 < q0) return ML_ERR_OVERFLOW;
        p0 = p1;
        p1 = p2;
        q0 = q1;
        q1 = q2;
        if (p1 != 0 && p1 > UINT64_MAX / p1) return ML_ERR_OVERFLOW;
        lhs = p1 * p1;
        if (q1 != 0 && D > UINT64_MAX / q1) return ML_ERR_OVERFLOW;
        {
            uint64_t dq = D * q1;
            if (q1 != 0 && dq > UINT64_MAX / q1) return ML_ERR_OVERFLOW;
            rhs = dq * q1;
        }
        if (rhs != UINT64_MAX && lhs == rhs + 1) {
            *x = p1;
            *y = q1;
            return ML_SUCCESS;
        }
        if (a > UINT64_MAX - m) return ML_ERR_OVERFLOW;
        m2 = d * a - m;
        if (m2 * m2 > D) return ML_ERR_OVERFLOW;
        d2 = (D - m2 * m2) / d;
        if (d2 == 0) return ML_ERR_INTERNAL;
        a = (a0 + m2) / d2;
        m = m2;
        d = d2;
    }
    return ML_ERR_UNSUPPORTED;
}
