#include "ml_compiler.h"
#include "ml_combinatorics.h"
#include <stdint.h>

/* v11S CLOSURE IP-12: portable combinatorics + overflow detection */

static uint64_t ml_gcd_u64(uint64_t a, uint64_t b) {
    while (b != 0) {
        uint64_t t = a % b;
        a = b;
        b = t;
    }
    return a;
}

ML_API uint64_t ml_gcd(uint64_t a, uint64_t b) {
    return ml_gcd_u64(a, b);
}

ML_API uint64_t ml_lcm(uint64_t a, uint64_t b) {
    if (a == 0 || b == 0) return 0;
    {
        uint64_t g = ml_gcd_u64(a, b);
        uint64_t d = a / g;
        if (d > UINT64_MAX / b) return UINT64_MAX;
        return d * b;
    }
}

ML_API uint64_t ml_factorial(int x) {
    if (x < 0) return UINT64_MAX;
    if (x > 20) return UINT64_MAX;
    if (x == 0) return 1;

    uint64_t result = 1;

    for (int i = 2; i <= x; i++) {
        uint64_t factor = (uint64_t)i;

        if (result > UINT64_MAX / factor) {
            return UINT64_MAX;
        }

        result *= factor;
    }

    return result;
}

ML_API uint64_t ml_npr(int n, int r) {
    if (n < 0 || r < 0 || r > n) return 0;

    uint64_t result = 1;

    for (int i = 0; i < r; i++) {
        uint64_t factor = (uint64_t)n - (uint64_t)i;

        if (factor == 0) return 0;

        if (result > UINT64_MAX / factor) {
            return UINT64_MAX;
        }

        result *= factor;
    }

    return result;
}

ML_API uint64_t ml_ncr(int n, int r) {
    if (n < 0 || r < 0 || r > n) return 0;

    if (r > n - r) {
        r = n - r;
    }

    /*
     * DoS guard: r = min(k,n-k) loop is O(r) with GCD per step.
     * For r > 1024 with n >= 2r, C(n,r) >= C(2r,r) >> UINT64_MAX
     * (C(2048,1024) ~ 10^613), so overflow is certain. Return the
     * overflow sentinel immediately instead of looping up to ~1B times
     * (e.g. ml_ncr(INT_MAX, INT_MAX/2)).
     */
    if (r > 1024) {
        return UINT64_MAX;
    }

    uint64_t result = 1;

    for (int i = 1; i <= r; i++) {
        uint64_t a = (uint64_t)n - (uint64_t)r + (uint64_t)i;
        uint64_t b = (uint64_t)i;

        /*
         * Cancel factors before multiplication.
         *
         * This keeps the intermediate result as small as possible while
         * preserving exact integer arithmetic.
         */
        uint64_t g = ml_gcd_u64(result, b);
        result /= g;
        b /= g;

        g = ml_gcd_u64(a, b);
        a /= g;
        b /= g;

        /*
         * For valid binomial coefficients, b should now be 1.
         *
         * This defensive path exists only to guarantee termination and
         * correctness if an unexpected cancellation remainder appears.
         */
        if (b != 1) {
            if (a % b != 0) {
                return UINT64_MAX;
            }
            a /= b;
        }

        if (a != 1) {
            if (result > UINT64_MAX / a) {
                return UINT64_MAX;
            }
            result *= a;
        }
    }

    return result;
}

static uint64_t ml_addmod(uint64_t a, uint64_t b, uint64_t m) {
    /* (a+b)%m for a,b<m without overflow: m-a never overflows. */
    uint64_t d = m - a;
    if (b >= d) return b - d;
    return a + b;
}

static uint64_t ml_mulmod(uint64_t a, uint64_t b, uint64_t m) {
    uint64_t r = 0;
    a %= m;
    while (b) {
        if (b & 1ULL) r = ml_addmod(r, a, m);
        a = ml_addmod(a, a, m);
        b >>= 1;
    }
    return r;
}

ML_API uint64_t ml_modpow(uint64_t a, uint64_t e, uint64_t m) {
    if (m <= 1) return 0;
    {
        uint64_t r = 1 % m;
        a %= m;
        while (e) {
            if (e & 1ULL) r = ml_mulmod(r, a, m);
            a = ml_mulmod(a, a, m);
            e >>= 1;
        }
        return r;
    }
}

ML_API int64_t ml_egcd(int64_t a, int64_t b, int64_t *x, int64_t *y) {
    /* INT64_MIN guards: -a overflows and MIN%-1 SIGFPEs. Use unsigned
     * magnitude path for the b==0 base case; reject MIN/-1 recursion. */
    if (a == INT64_MIN || b == INT64_MIN) {
        if (x) *x = 0;
        if (y) *y = 0;
        return -1;
    }
    if (b == 0) {
        if (x) *x = (a >= 0) ? 1 : -1;
        if (y) *y = 0;
        return (a >= 0) ? a : -a;
    }
    {
        int64_t x1 = 0, y1 = 0;
        int64_t g = ml_egcd(b, a % b, &x1, &y1);
        if (x) *x = y1;
        if (y) *y = x1 - (a / b) * y1;
        return g;
    }
}

ML_API uint64_t ml_modinv(uint64_t a, uint64_t m) {
    if (m <= 1) return UINT64_MAX;
    /* int64 Euclid only covers m<=INT64_MAX; larger moduli would wrap
     * negative and give wrong residues. */
    if (m > (uint64_t)INT64_MAX) return UINT64_MAX;
    if (a >= m) a %= m;
    {
        int64_t x = 0, y = 0;
        int64_t g = ml_egcd((int64_t)(a % m), (int64_t)m, &x, &y);
        if (g != 1 && g != -1) return UINT64_MAX;
        {
            int64_t r = x % (int64_t)m;
            if (r < 0) r += (int64_t)m;
            return (uint64_t)r;
        }
    }
}

static int ml_mr_pass(uint64_t n, uint64_t a) {
    uint64_t d = n - 1;
    int s = 0;
    while ((d & 1ULL) == 0ULL) { d >>= 1; s++; }
    {
        uint64_t x = ml_modpow(a % n, d, n);
        if (x == 1 || x == n - 1) return 1;
        for (int r = 1; r < s; r++) {
            x = ml_mulmod(x, x, n);
            if (x == n - 1) return 1;
        }
        return 0;
    }
}

ML_API int ml_is_prime64(uint64_t n) {
    static const uint64_t small[] = {2,3,5,7,11,13,17,19,23,29,31,37};
    if (n < 2) return 0;
    for (unsigned i = 0; i < sizeof(small)/sizeof(small[0]); i++) {
        if (n == small[i]) return 1;
        if (n % small[i] == 0ULL) return 0;
    }
    /* Deterministic MR for 64-bit (Sinclair bases). */
    static const uint64_t bases[] = {2ULL, 325ULL, 9375ULL, 28178ULL, 450775ULL, 9780504ULL, 1795265022ULL};
    for (unsigned i = 0; i < sizeof(bases)/sizeof(bases[0]); i++) {
        uint64_t a = bases[i] % n;
        if (a == 0) continue;
        if (!ml_mr_pass(n, a)) return 0;
    }
    return 1;
}

static uint64_t ml_rho_f(uint64_t x, uint64_t c, uint64_t m) {
    return (ml_mulmod(x, x, m) + c) % m;
}

static uint64_t ml_pollard_rho(uint64_t n) {
    if (n % 2ULL == 0ULL) return 2ULL;
    {
        uint64_t x = 2ULL, y = 2ULL, d = 1ULL, c = 1ULL;
        while (d == 1ULL) {
            x = ml_rho_f(x, c, n);
            y = ml_rho_f(ml_rho_f(y, c, n), c, n);
            {
                uint64_t diff = (x > y) ? x - y : y - x;
                if (diff == 0) { c++; x = y = 2ULL; continue; }
                /* gcd via Euclid (64-bit, no recursion). */
                uint64_t a = diff, b = n;
                while (b) { uint64_t t = a % b; a = b; b = t; }
                d = a;
                if (d == n) { c++; x = y = 2ULL; d = 1ULL; }
            }
        }
        return d;
    }
}

static void ml_factor_rec(uint64_t n, uint64_t *out, int *cnt, int cap) {
    if (n <= 1 || *cnt >= cap) return;
    if (ml_is_prime64(n)) {
        /* Insert sorted. */
        int i = *cnt;
        while (i > 0 && out[i - 1] > n) { out[i] = out[i - 1]; i--; }
        out[i] = n;
        (*cnt)++;
        return;
    }
    {
        uint64_t d = ml_pollard_rho(n);
        if (d <= 1 || d >= n) {
            /* Fallback trial division for small factors. */
            for (uint64_t p = 3; p * p <= n && p < 1000000ULL; p += 2) {
                if (n % p == 0ULL) { d = p; break; }
            }
            if (d <= 1 || d >= n) return;
        }
        ml_factor_rec(d, out, cnt, cap);
        ml_factor_rec(n / d, out, cnt, cap);
    }
}

ML_API uint64_t ml_phi(uint64_t n) {
    if (n == 0) return 0;
    if (n == 1) return 1;
    {
        uint64_t fac[64];
        int cnt = 0;
        ml_factor_rec(n, fac, &cnt, 64);
        if (cnt == 0) return 0;
        {
            uint64_t r = n;
            uint64_t prev = 0;
            for (int i = 0; i < cnt; i++) {
                if (fac[i] != prev) {
                    r = r / fac[i] * (fac[i] - 1ULL);
                    prev = fac[i];
                }
            }
            return r;
        }
    }
}

ML_API uint64_t ml_catalan(int n) {
    if (n < 0) return 0;
    if (n == 0) return 1;
    /* C37 first overflows 64-bit; C35/C36 fit but C(70,35)/C(72,36) do not,
     * so ml_ncr(2n,n)/(n+1) falsely overflows. DESPOT-FIX: recurrence
     * C_{k+1}=C_k*2(2k+1)/(k+2) with GCD cancellation keeps every
     * intermediate <= 1.2e19 (C99, no __int128). */
    if (n > 36) return UINT64_MAX;
    {
        uint64_t c = 1;
        for (int k = 0; k < n; k++) {
            uint64_t num = 2u * ((uint64_t)k * 2u + 1u);
            /* num = 2*(2k+1); den = k+2 */
            uint64_t den = (uint64_t)k + 2u;
            uint64_t g = ml_gcd_u64(c, den);
            c /= g;
            den /= g;
            g = ml_gcd_u64(num, den);
            num /= g;
            den /= g;
            if (den != 1) {
                if (num % den != 0) return UINT64_MAX;
                num /= den;
            }
            if (num != 1 && c > UINT64_MAX / num) return UINT64_MAX;
            c *= num;
        }
        return c;
    }
}


ML_API uint64_t ml_fib_pair(long long n, uint64_t *fn, uint64_t *fn1) {
    if (n < 0) return 0;
    /* F(93) fits, F(94) overflows; O(n) loop would hang on 1e12. */
    if (n > 93) {
        if (fn) *fn = UINT64_MAX;
        if (fn1) *fn1 = UINT64_MAX;
        return UINT64_MAX;
    }
    /* Iterative O(n), n<=93 fits; overflow -> MAX sentinel. */
    if (n == 0) {
        if (fn) *fn = 0;
        if (fn1) *fn1 = 1;
        return 0;
    }
    {
        uint64_t a = 0, b = 1;
        for (long long i = 1; i < n; i++) {
            if (b > UINT64_MAX - a) {
                if (fn) *fn = UINT64_MAX;
                if (fn1) *fn1 = UINT64_MAX;
                return UINT64_MAX;
            }
            {
                uint64_t c = a + b;
                a = b; b = c;
            }
        }
        if (fn) *fn = b;
        if (fn1) {
            if (b > UINT64_MAX - a) *fn1 = UINT64_MAX;
            else *fn1 = a + b;
        }
        return b;
    }
}

ML_API uint64_t ml_derange(int n) {
    if (n < 0) return 0;
    if (n == 0) return 1;
    if (n == 1) return 0;
    if (n > 20) return UINT64_MAX;
    {
        uint64_t a = 1, b = 0, c = 0;
        for (int i = 2; i <= n; i++) {
            uint64_t s = a + b;
            uint64_t m;
            if (s < a) return UINT64_MAX;
            m = (uint64_t)(i - 1);
            if (m != 0 && s > UINT64_MAX / m) return UINT64_MAX;
            c = s * m;
            /* exact overflow check for s*(i-1): */
            if ((uint64_t)(i - 1) != 0 && c / (uint64_t)(i - 1) != s) return UINT64_MAX;
            a = b; b = c;
        }
        return b;
    }
}
