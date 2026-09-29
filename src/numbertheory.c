#include "ml_compiler.h"
#include "ml_numbertheory.h"
#include "ml_combinatorics.h"
#include <stdlib.h>

static uint64_t ml_abs_u64(int64_t v) {
    if (v >= 0) return (uint64_t)v;
    if (v == INT64_MIN) return (uint64_t)INT64_MAX + 1ULL;
    return (uint64_t)(-v);
}

ML_API int ml_legendre_symbol(int64_t a, uint64_t p) {
    if (p < 3 || (p & 1ULL) == 0ULL) return 0;
    if (!ml_is_prime64(p)) return 0;
    {
        int64_t pp = (int64_t)p;
        int64_t t = a % pp;
        if (t < 0) t += pp;
        uint64_t am = (uint64_t)t;
        if (am == 0) return 0;
        uint64_t r = ml_modpow(am, (p - 1) / 2, p);
        if (r == 1) return 1;
        if (r == p - 1) return -1;
        return 0;
    }
}

ML_API int ml_jacobi_symbol(int64_t a, uint64_t n) {
    if (n == 0 || (n & 1ULL) == 0ULL) return 0;
    if (n == 1) return 1;
    {
        uint64_t aa, nn = n;
        int t = 1;
        /* (a/n) depends only on a mod n, reduced into [0,n).  Reducing in
         * the magnitude domain keeps this exact for a = INT64_MIN and for
         * n > 2^63-1, and the loop below stays entirely unsigned. */
        if (a < 0) {
            uint64_t m = (uint64_t)(-(a + 1)) + 1ULL;   /* |a|, no INT64_MIN UB */
            uint64_t r = m % nn;
            aa = (r == 0ULL) ? 0ULL : nn - r;
        } else {
            aa = (uint64_t)a % nn;
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

ML_API int ml_kronecker_symbol(int64_t a, int64_t n) {
    if (n == 0) {
        uint64_t aa = ml_abs_u64(a);
        return (aa == 1) ? 1 : 0;
    }
    if (n == 1) return 1;
    if (n == -1) return (a < 0) ? -1 : 1;
    {
        int t = 1;
        int64_t nn = n;
        if (nn < 0) {
            nn = -nn;
            if (a < 0) t = -t;
        }
        /* Factor out powers of 2: (a/2) = 0 if even, else +-1 by a mod 8. */
        int e = 0;
        while ((nn & 1LL) == 0) { nn >>= 1; e++; }
        if (e > 0) {
            uint64_t am8 = ml_abs_u64(a) & 7ULL;
            if ((a & 1LL) == 0) return 0;
            int a2 = 0;
            if (am8 == 1 || am8 == 7) a2 = 1;
            else if (am8 == 3 || am8 == 5) a2 = -1;
            if (e % 2 == 1) t *= a2;
            else if (a2 == 0) return 0;
        }
        if (nn == 1) return t;
        return t * ml_jacobi_symbol(a, (uint64_t)nn);
    }
}

ML_API int ml_mobius(int64_t n) {
    if (n == 0) return 0;
    uint64_t m = ml_abs_u64(n);
    if (m == 1) return 1;
    {
        int prime_factors = 0;
        uint64_t tmp = m;
        for (uint64_t p = 2; p <= tmp / p; p += (p == 2 ? 1 : 2)) {
            if (p > 1000000ULL) {
                /* Too large for trial: use prime test to decide. If the
                 * cofactor is prime, count it; else unknown -> 2. */
                if (ml_is_prime64(tmp)) { prime_factors++; tmp = 1; break; }
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

ML_API uint64_t ml_mult_order(uint64_t a, uint64_t m) {
    if (m <= 1) return 0;
    a %= m;
    if (ml_gcd(a, m) != 1) return 0;
    {
        uint64_t phi = ml_phi(m);
        if (phi == 0) return 0;
        /* Order divides phi: strip prime factors (trial to 1e6, else phi). */
        uint64_t ord = phi, tmp = phi;
        for (uint64_t p = 2; p <= tmp / p; p += (p == 2 ? 1 : 2)) {
            if (p > 1000000ULL) break;
            if (tmp % p == 0) {
                while (tmp % p == 0) tmp /= p;
                while (ord % p == 0 && ml_modpow(a, ord / p, m) == 1) ord /= p;
            }
        }
        if (tmp > 1 && ord % tmp == 0 && ml_modpow(a, ord / tmp, m) == 1) ord /= tmp;
        return ord;
    }
}

ML_API int ml_is_primitive_root(uint64_t a, uint64_t p) {
    if (p < 2 || !ml_is_prime64(p)) return 0;
    a %= p;
    if (a == 0) return 0;
    return ml_mult_order(a, p) == p - 1 ? 1 : 0;
}

ML_API ml_status_t ml_crt2(uint64_t a1, uint64_t m1, uint64_t a2, uint64_t m2,
                           uint64_t *x, uint64_t *lcm) {
    if (ML_UNLIKELY(!x || !lcm || m1 == 0 || m2 == 0)) return ML_ERR_INVALID_ARG;
    {
        uint64_t g = ml_gcd(m1, m2);
        uint64_t d1 = (a2 >= a1) ? a2 - a1 : a1 - a2;
        if (d1 % g != 0) return ML_ERR_SINGULAR;
        {
            uint64_t l = m1 / g * m2;
            if (l / m2 != m1 / g) return ML_ERR_WORKSPACE;
            /* Solve m1*t = a2-a1 (mod m2): reduce by g, invert. */
            uint64_t p = m1 / g, q = m2 / g;
            uint64_t diff = (a2 >= a1) ? (a2 - a1) / g : (m2 - ((a1 - a2) / g) % m2) % m2;
            (void)q;
            {
                uint64_t inv = ml_modinv(p % m2, m2 / g == 0 ? 1 : q);
                if (inv == UINT64_MAX) {
                    /* Fallback extended Euclid on reduced moduli. */
                    int64_t ex = 0, ey = 0;
                    int64_t gg = ml_egcd((int64_t)p, (int64_t)q, &ex, &ey);
                    if (gg != 1 && gg != -1) return ML_ERR_SINGULAR;
                    ex %= (int64_t)q;
                    if (ex < 0) ex += (int64_t)q;
                    inv = (uint64_t)ex;
                }
                __uint128_t t = (__uint128_t)diff * inv % q;
                __uint128_t r = (__uint128_t)a1 + (__uint128_t)m1 * t;
                *x = (uint64_t)(r % l);
                *lcm = l;
                return ML_SUCCESS;
            }
        }
    }
}

ML_API uint64_t ml_prime_pi(uint64_t n) {
    if (n < 2) return 0;
    if (n > 100000000ULL) return UINT64_MAX;
    {
        uint64_t nn = n;
        /* Heap-backed odd-only sieve sized by the limit: avoids a 50MB
         * static BSS image and is thread-safe (no shared mutable state). */
        uint64_t size = nn / 2;
        uint8_t *bits = (uint8_t *)malloc((size_t)size * sizeof(uint8_t));
        if (!bits) return UINT64_MAX;
        for (uint64_t i = 0; i < size; i++) bits[i] = 1;
        for (uint64_t p = 3; p * p <= nn; p += 2) {
            if (bits[p / 2]) {
                for (uint64_t q2 = p * p; q2 <= nn; q2 += 2 * p) bits[q2 / 2] = 0;
            }
        }
        uint64_t c = 1;
        for (uint64_t i = 1; i < size; i++) {
            if (bits[i]) c++;
        }
        free(bits);
        return c;
    }
}

ML_API int ml_quadratic_reciprocity_check(uint64_t p, uint64_t q) {
    if (!ml_is_prime64(p) || !ml_is_prime64(q) || p < 3 || q < 3 || p == q) return 0;
    {
        int64_t lp = ml_legendre_symbol((int64_t)p, q);
        int64_t lq = ml_legendre_symbol((int64_t)q, p);
        int sign = (((p - 1) / 2) & ((q - 1) / 2) & 1) ? -1 : 1;
        return (lp == sign * lq) ? 1 : 0;
    }
}
