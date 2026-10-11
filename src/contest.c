#include "ml_compiler.h"
#include "ml_contest.h"
#include "ml_core.h"
#include "ml_combinatorics.h"
#include "ml_integral.h"
#include <stddef.h>
#include <stdlib.h>

/* Contest mathematics: recreational -> olympiad -> Putnam.
 * Exact u64 paths use overflow-checked arithmetic with UINT64_MAX sentinel
 * (house convention); double paths are exactness-guarded (sum_k pattern)
 * or documented LD approximations (regression-guard tier). */

static uint64_t mlc_gcd_u64(uint64_t a, uint64_t b) {
    while (b != 0) { uint64_t t = a % b; a = b; b = t; }
    return a;
}

static uint64_t mlc_isqrt_u64(uint64_t x) {
    /* Integer Newton: no FP involved, exact for the full u64 range. */
    uint64_t r;
    uint64_t y;
    if (x <= 1) return x;
    r = x;
    y = (r + 1) / 2;
    while (y < r) { r = y; y = (r + x / r) / 2; }
    return r;
}

static uint64_t mlc_addmod_u64(uint64_t a, uint64_t b, uint64_t m) {
    uint64_t d = m - a;
    if (b >= d) return b - d;
    return a + b;
}

static uint64_t mlc_mulmod_u64(uint64_t a, uint64_t b, uint64_t m) {
    /* Portable double-and-add: no __int128 (strict-C99 pedantic). */
    uint64_t r = 0;
    a %= m;
    while (b != 0) {
        if ((b & 1ULL) != 0) r = mlc_addmod_u64(r, a, m);
        a = mlc_addmod_u64(a, a, m);
        b >>= 1;
    }
    return r;
}

/* ---- Tier 1 ---- */

ML_API uint64_t ml_triangular(uint64_t n) {
    /* T_n = n(n+1)/2: halve first, then checked multiply. */
    uint64_t a = n, b = n + 1;
    if (b < a) return UINT64_MAX;
    if ((a & 1ULL) == 0) a /= 2;
    else b /= 2;
    if (a != 0 && b > UINT64_MAX / a) return UINT64_MAX;
    return a * b;
}

ML_API uint64_t ml_pentagonal(uint64_t n) {
    /* P_n = n(3n-1)/2. */
    uint64_t t;
    uint64_t a;
    uint64_t b;
    if (n != 0 && 3 > UINT64_MAX / n) return UINT64_MAX;
    t = 3 * n;
    if (t < 1) return UINT64_MAX;
    a = n;
    b = t - 1;
    if ((a & 1ULL) == 0) a /= 2;
    else b /= 2;
    if (a != 0 && b > UINT64_MAX / a) return UINT64_MAX;
    return a * b;
}

ML_API uint64_t ml_hexagonal(uint64_t n) {
    /* H_n = n(2n-1). */
    uint64_t b;
    if (n != 0 && 2 > UINT64_MAX / n) return UINT64_MAX;
    b = 2 * n;
    if (b < 1) return UINT64_MAX;
    b -= 1;
    if (n != 0 && b > UINT64_MAX / n) return UINT64_MAX;
    return n * b;
}

ML_API int ml_is_triangular(uint64_t m) {
    /* m triangular iff 8m+1 is an odd square. */
    uint64_t t;
    uint64_t s;
    if (m > (UINT64_MAX - 1ULL) / 8ULL) return 0;
    t = 8ULL * m + 1ULL;
    s = mlc_isqrt_u64(t);
    return (s * s == t) ? 1 : 0;
}

ML_API uint64_t ml_double_factorial(long long n) {
    /* n!! with (-1)!! = 0!! = 1 convention. */
    uint64_t r = 1;
    long long k;
    if (n < -1) return UINT64_MAX;
    if (n <= 1) return 1;
    for (k = n; k >= 2; k -= 2) {
        uint64_t f = (uint64_t)k;
        if (r > UINT64_MAX / f) return UINT64_MAX;
        r *= f;
    }
    return r;
}

ML_API int ml_binom_odd(uint64_t n, uint64_t k) {
    /* Lucas mod 2: C(n,k) odd iff every bit of k is in n. */
    if (k > n) return 0;
    return ((k & ~n) == 0) ? 1 : 0;
}

ML_API uint64_t ml_stars_bars(uint64_t stars, uint64_t bars) {
    /* C(stars+bars, bars); ml_ncr is int-based, so cap at INT_MAX. */
    uint64_t tot;
    if (stars > UINT64_MAX - bars) return UINT64_MAX;
    tot = stars + bars;
    if (tot > (uint64_t)2147483647) return UINT64_MAX;
    if (bars > (uint64_t)2147483647) return UINT64_MAX;
    return ml_ncr((int)tot, (int)bars);
}

ML_API uint64_t ml_hanoi_moves(int n) {
    /* 2^n - 1; n = 64 would equal the sentinel, so it is capped. */
    if (n < 0 || n >= 64) return UINT64_MAX;
    if (n == 0) return 0;
    return (1ULL << (unsigned)n) - 1ULL;
}

ML_API uint64_t ml_josephus(uint64_t n, uint64_t k) {
    /* 0-indexed survivor, iterative. O(n): cap keeps worst case ~0.3s. */
    uint64_t j = 0;
    uint64_t i;
    if (n < 1 || k < 1) return UINT64_MAX;
    if (n > 100000000ULL) return UINT64_MAX;
    for (i = 2; i <= n; i++) j = (j + k) % i;
    return j;
}

ML_API uint64_t ml_lucas_number(int n) {
    /* L_0=2, L_1=1. L_93 fits, L_94 overflows: checked loop decides. */
    uint64_t a = 2, b = 1;
    int i;
    if (n < 0 || n > 1000000) return UINT64_MAX;
    if (n == 0) return 2;
    if (n == 1) return 1;
    for (i = 2; i <= n; i++) {
        uint64_t c;
        if (b > UINT64_MAX - a) return UINT64_MAX;
        c = a + b;
        a = b;
        b = c;
    }
    return b;
}

ML_API uint64_t ml_pell_number(int n) {
    /* P_0=0, P_1=1, P_n=2P_{n-1}+P_{n-2}. */
    uint64_t a = 0, b = 1;
    int i;
    if (n < 0 || n > 1000000) return UINT64_MAX;
    if (n == 0) return 0;
    if (n == 1) return 1;
    for (i = 2; i <= n; i++) {
        uint64_t c;
        uint64_t t;
        if (b > (UINT64_MAX - a) / 2) return UINT64_MAX;
        t = 2 * b;
        if (t > UINT64_MAX - a) return UINT64_MAX;
        c = t + a;
        a = b;
        b = c;
    }
    return b;
}

ML_API uint64_t ml_tribonacci(int n) {
    /* T_0=0, T_1=0, T_2=1. */
    uint64_t a = 0, b = 0, c = 1;
    int i;
    if (n < 0 || n > 1000000) return UINT64_MAX;
    if (n == 0 || n == 1) return 0;
    if (n == 2) return 1;
    for (i = 3; i <= n; i++) {
        uint64_t d;
        if (a > UINT64_MAX - b) return UINT64_MAX;
        d = a + b;
        if (d > UINT64_MAX - c) return UINT64_MAX;
        d += c;
        a = b;
        b = c;
        c = d;
    }
    return c;
}

/* ---- Tier 2 ---- */

ML_API uint64_t ml_bell_number(int n) {
    /* Aitken's array: row[0] = previous last; row[j] = row[j-1]+old[j-1].
     * B_25 fits, B_26 overflows. */
    uint64_t row[26];
    uint64_t last = 1;
    int i;
    int j;
    if (n < 0 || n > 25) return (n == 0) ? 1 : UINT64_MAX;
    if (n == 0) return 1;
    for (i = 0; i < 26; i++) row[i] = 0;
    row[0] = 1;
    for (i = 1; i <= n; i++) {
        uint64_t prev = row[0];
        row[0] = last;
        for (j = 1; j <= i; j++) {
            uint64_t old = row[j];
            if (row[j - 1] > UINT64_MAX - prev) return UINT64_MAX;
            row[j] = row[j - 1] + prev;
            prev = old;
        }
        last = row[i];
    }
    return row[0];
}

ML_API uint64_t ml_stirling2(int n, int k) {
    /* S(n,k) = k*S(n-1,k) + S(n-1,k-1). Table 31x31 fits stack. */
    uint64_t s[31][31];
    int i;
    int j;
    if (n < 0 || k < 0 || k > n || n > 30) return (n == 0 && k == 0) ? 1 : UINT64_MAX;
    for (i = 0; i <= n; i++)
        for (j = 0; j <= n; j++) s[i][j] = 0;
    s[0][0] = 1;
    for (i = 1; i <= n; i++) {
        for (j = 1; j <= i; j++) {
            uint64_t t1 = 0, t2;
            uint64_t kj = (uint64_t)j;
            if (s[i - 1][j] != 0) {
                if (kj != 0 && s[i - 1][j] > UINT64_MAX / kj) return UINT64_MAX;
                t1 = kj * s[i - 1][j];
            }
            t2 = s[i - 1][j - 1];
            if (t1 > UINT64_MAX - t2) return UINT64_MAX;
            s[i][j] = t1 + t2;
        }
    }
    return s[n][k];
}

ML_API uint64_t ml_stirling1(int n, int k) {
    /* Unsigned |s(n,k)| = |s(n-1,k-1)| + (n-1)*|s(n-1,k)|. */
    uint64_t s[31][31];
    int i;
    int j;
    if (n < 0 || k < 0 || k > n || n > 30) return (n == 0 && k == 0) ? 1 : UINT64_MAX;
    for (i = 0; i <= n; i++)
        for (j = 0; j <= n; j++) s[i][j] = 0;
    s[0][0] = 1;
    for (i = 1; i <= n; i++) {
        for (j = 1; j <= i; j++) {
            uint64_t t1 = 0, t2;
            uint64_t nm1 = (uint64_t)(i - 1);
            if (s[i - 1][j] != 0) {
                if (nm1 != 0 && s[i - 1][j] > UINT64_MAX / nm1) return UINT64_MAX;
                t1 = nm1 * s[i - 1][j];
            }
            t2 = s[i - 1][j - 1];
            if (t1 > UINT64_MAX - t2) return UINT64_MAX;
            s[i][j] = t1 + t2;
        }
    }
    return s[n][k];
}

ML_API uint64_t ml_eulerian(int n, int m) {
    /* A(n,m) = (n-m)*A(n-1,m-1) + (m+1)*A(n-1,m), 0 <= m < n. */
    uint64_t a[31][31];
    int i;
    int j;
    if (n < 0 || m < 0 || n > 30) return UINT64_MAX;
    if (n == 0) return (m == 0) ? 1 : 0;
    if (m >= n) return 0;
    for (i = 0; i <= n; i++)
        for (j = 0; j <= n; j++) a[i][j] = 0;
    a[0][0] = 1;
    for (i = 1; i <= n; i++) {
        for (j = 0; j < i; j++) {
            uint64_t t1 = 0, t2 = 0;
            uint64_t c1 = (uint64_t)(i - j);
            uint64_t c2 = (uint64_t)(j + 1);
            if (j > 0 && a[i - 1][j - 1] != 0) {
                if (c1 != 0 && a[i - 1][j - 1] > UINT64_MAX / c1) return UINT64_MAX;
                t1 = c1 * a[i - 1][j - 1];
            }
            if (a[i - 1][j] != 0) {
                if (c2 != 0 && a[i - 1][j] > UINT64_MAX / c2) return UINT64_MAX;
                t2 = c2 * a[i - 1][j];
            }
            if (t1 > UINT64_MAX - t2) return UINT64_MAX;
            a[i][j] = t1 + t2;
        }
    }
    return a[n][m];
}

ML_API uint64_t ml_lah_number(int n, int k) {
    /* L(n,k) = C(n-1,k-1) * n! / k!. */
    uint64_t c;
    uint64_t fn;
    uint64_t fk;
    uint64_t r;
    if (n <= 0 || k <= 0 || k > n) return 0;
    if (n - 1 > 2147483647 || k - 1 > 2147483647) return UINT64_MAX;
    c = ml_ncr(n - 1, k - 1);
    fn = ml_factorial(n);
    fk = ml_factorial(k);
    if (c == UINT64_MAX || fn == UINT64_MAX || fk == UINT64_MAX) return UINT64_MAX;
    if (c != 0 && fn > UINT64_MAX / c) return UINT64_MAX;
    r = c * fn;
    if (fk == 0 || r % fk != 0) return UINT64_MAX;
    return r / fk;
}

ML_API uint64_t ml_motzkin(int n) {
    /* (n+2)M_n = (2n+1)M_{n-1} + 3(n-1)M_{n-2}; exact division. */
    uint64_t m0 = 1, m1 = 1;
    int i;
    if (n < 0 || n > 35) return (n == 0) ? 1 : UINT64_MAX;
    if (n == 0 || n == 1) return 1;
    for (i = 2; i <= n; i++) {
        uint64_t u = (uint64_t)(2 * i + 1);
        uint64_t v = (uint64_t)(3 * (i - 1));
        uint64_t den = (uint64_t)(i + 2);
        uint64_t t1;
        uint64_t t2;
        uint64_t num;
        if (m1 != 0 && u > UINT64_MAX / m1) return UINT64_MAX;
        t1 = u * m1;
        if (m0 != 0 && v > UINT64_MAX / m0) return UINT64_MAX;
        t2 = v * m0;
        if (t1 > UINT64_MAX - t2) return UINT64_MAX;
        num = t1 + t2;
        if (den == 0 || num % den != 0) return UINT64_MAX;
        m0 = m1;
        m1 = num / den;
    }
    return m1;
}

ML_API uint64_t ml_schroder_large(int n) {
    /* (n+1)S_n = 3(2n-1)S_{n-1} - (n-2)S_{n-2}. */
    uint64_t s0 = 1, s1 = 2;
    int i;
    if (n < 0 || n > 25) return (n == 0) ? 1 : UINT64_MAX;
    if (n == 0) return 1;
    if (n == 1) return 2;
    for (i = 2; i <= n; i++) {
        uint64_t u = (uint64_t)(3 * (2 * i - 1));
        uint64_t v = (uint64_t)(i - 2);
        uint64_t den = (uint64_t)(i + 1);
        uint64_t t1;
        uint64_t t2;
        uint64_t num;
        if (s1 != 0 && u > UINT64_MAX / s1) return UINT64_MAX;
        t1 = u * s1;
        if (v != 0 && s0 > UINT64_MAX / v) return UINT64_MAX;
        t2 = v * s0;
        if (t1 < t2) return UINT64_MAX;
        num = t1 - t2;
        if (den == 0 || num % den != 0) return UINT64_MAX;
        s0 = s1;
        s1 = num / den;
    }
    return s1;
}

ML_API uint64_t ml_narayana(int n, int k) {
    /* N(n,k) = C(n,k)C(n,k-1)/n, 1 <= k <= n. Sum_k N = Catalan_n. */
    uint64_t c1;
    uint64_t c2;
    uint64_t r;
    uint64_t un;
    if (n <= 0 || k <= 0 || k > n) return 0;
    if (n > 67) return UINT64_MAX;
    c1 = ml_ncr(n, k);
    c2 = ml_ncr(n, k - 1);
    if (c1 == UINT64_MAX || c2 == UINT64_MAX) return UINT64_MAX;
    if (c1 != 0 && c2 > UINT64_MAX / c1) return UINT64_MAX;
    r = c1 * c2;
    un = (uint64_t)n;
    if (un == 0 || r % un != 0) return UINT64_MAX;
    return r / un;
}

ML_API uint64_t ml_part_distinct(int n) {
    /* Coin DP over distinct parts: dp[i] += dp[i-k]. Euler: == p_odd(n). */
    uint64_t dp[201];
    int kk;
    int ii;
    if (n < 0 || n > 200) return (n == 0) ? 1 : UINT64_MAX;
    for (ii = 0; ii <= 200; ii++) dp[ii] = 0;
    dp[0] = 1;
    for (kk = 1; kk <= n; kk++) {
        for (ii = n; ii >= kk; ii--) {
            if (dp[ii - kk] > UINT64_MAX - dp[ii]) return UINT64_MAX;
            dp[ii] += dp[ii - kk];
        }
    }
    return dp[n];
}

ML_API uint64_t ml_goldbach_partitions(uint64_t n) {
    /* Unordered prime pairs p <= q, p+q = n. Heap odd sieve to n. */
    uint64_t size;
    uint8_t *bits;
    uint64_t p;
    uint64_t cnt = 0;
    uint64_t q;
    if (n < 4 || (n & 1ULL) != 0) return 0;
    if (n > 1000000ULL) return UINT64_MAX;
    size = (n + 1) / 2;
    if ((size_t)size != size) return UINT64_MAX;
    bits = (uint8_t *)malloc((size_t)size * sizeof(uint8_t));
    if (!bits) return UINT64_MAX;
    for (q = 0; q < size; q++) bits[q] = 1;
    for (p = 3; p * p <= n; p += 2) {
        if (bits[p / 2]) {
            uint64_t qq;
            for (qq = p * p; qq <= n; qq += 2 * p) bits[(size_t)(qq / 2)] = 0;
        }
    }
    for (p = 2; p <= n / 2; p++) {
        uint64_t r = n - p;
        int pp = (p == 2) ? 1 : ((p & 1ULL) ? (int)bits[(size_t)(p / 2)] : 0);
        int rp = (r == 2) ? 1 : ((r & 1ULL) ? (int)bits[(size_t)(r / 2)] : 0);
        if (pp && rp) cnt++;
    }
    free(bits);
    return cnt;
}

ML_API int ml_is_sophie_germain(uint64_t p) {
    /* p and 2p+1 both prime; guard 2p+1 overflow. */
    if (p < 2) return 0;
    if (!ml_is_prime64(p)) return 0;
    if (p > (UINT64_MAX - 1ULL) / 2ULL) return 0;
    return ml_is_prime64(2 * p + 1) ? 1 : 0;
}

ML_API int ml_lucas_lehmer(int p) {
    /* Mersenne 2^p-1 prime test. 1 = prime, 0 = composite, -1 = invalid. */
    uint64_t m;
    uint64_t s;
    int i;
    if (p < 3 || p > 61) return -1;
    if (!ml_is_prime64((uint64_t)p)) return -1;
    m = (p == 64) ? UINT64_MAX : ((1ULL << (unsigned)p) - 1ULL);
    if (m < 3) return -1;
    s = 4 % m;
    for (i = 3; i <= p; i++) {
        uint64_t sq = mlc_mulmod_u64(s, s, m);
        s = (sq >= 2) ? sq - 2 : sq + m - 2;
    }
    return (s == 0) ? 1 : 0;
}

ML_API uint64_t ml_modsqrt_tonelli(uint64_t a, uint64_t p) {
    /* Least r with r^2 = a (mod p); UINT64_MAX when none/invalid. */
    uint64_t q;
    uint64_t s = 0;
    uint64_t z;
    uint64_t c;
    uint64_t r;
    uint64_t t;
    uint64_t m;
    if (p < 3 || (p & 1ULL) == 0) return UINT64_MAX;
    if (!ml_is_prime64(p)) return UINT64_MAX;
    a %= p;
    if (a == 0) return 0;
    if (ml_modpow(a, (p - 1) / 2, p) != 1) return UINT64_MAX;
    if ((p & 3ULL) == 3ULL) {
        r = ml_modpow(a, (p + 1) / 4, p);
        return (r <= p - r) ? r : p - r;
    }
    q = p - 1;
    while ((q & 1ULL) == 0) { q >>= 1; s++; }
    z = 2;
    while (ml_modpow(z, (p - 1) / 2, p) != p - 1) {
        z++;
        if (z >= p) return UINT64_MAX;
    }
    c = ml_modpow(z, q, p);
    r = ml_modpow(a, (q + 1) / 2, p);
    t = ml_modpow(a, q, p);
    m = s;
    while (t != 1) {
        uint64_t tt = t;
        uint64_t i;
        uint64_t b;
        uint64_t i_found = 0;
        for (i = 1; i < m; i++) {
            tt = mlc_mulmod_u64(tt, tt, p);
            if (tt == 1) { i_found = i; break; }
        }
        if (i_found == 0) return UINT64_MAX;
        b = c;
        for (i = 0; i < m - i_found - 1; i++) b = mlc_mulmod_u64(b, b, p);
        r = mlc_mulmod_u64(r, b, p);
        c = mlc_mulmod_u64(b, b, p);
        t = mlc_mulmod_u64(t, c, p);
        m = i_found;
    }
    return (r <= p - r) ? r : p - r;
}

ML_API uint64_t ml_tau_divisors(uint64_t n) {
    /* d(n) by trial to sqrt; cap 1e12 keeps worst loop at 1e6. */
    uint64_t cnt = 0;
    uint64_t d;
    if (n == 0) return 0;
    if (n > 1000000000000ULL) return UINT64_MAX;
    for (d = 1; d <= n / d; d++) {
        if (n % d == 0) cnt += (d * d == n) ? 1 : 2;
    }
    return cnt;
}

ML_API uint64_t ml_sigma_sum(uint64_t n) {
    /* sigma(n) with checked adds/muls. */
    uint64_t s = 0;
    uint64_t d;
    if (n == 0) return 0;
    if (n > 1000000000000ULL) return UINT64_MAX;
    for (d = 1; d <= n / d; d++) {
        if (n % d == 0) {
            uint64_t e = n / d;
            if (s > UINT64_MAX - d) return UINT64_MAX;
            s += d;
            if (e != d) {
                if (s > UINT64_MAX - e) return UINT64_MAX;
                s += e;
            }
        }
    }
    return s;
}

ML_API uint64_t ml_frobenius2(uint64_t a, uint64_t b) {
    /* Largest nonrepresentable by coprime a,b: ab-a-b. */
    uint64_t p;
    if (a < 2 || b < 2) return 0;
    if (mlc_gcd_u64(a, b) != 1) return 0;
    if (a != 0 && b > UINT64_MAX / a) return UINT64_MAX;
    p = a * b;
    if (p < a + b) return UINT64_MAX;
    return p - a - b;
}

ML_API uint64_t ml_frobenius_count2(uint64_t a, uint64_t b) {
    /* Count of nonrepresentable: (a-1)(b-1)/2. */
    uint64_t u;
    uint64_t v;
    if (a < 2 || b < 2) return 0;
    if (mlc_gcd_u64(a, b) != 1) return 0;
    u = a - 1;
    v = b - 1;
    if ((u & 1ULL) == 0) u /= 2;
    else v /= 2;
    if (u != 0 && v > UINT64_MAX / u) return UINT64_MAX;
    return u * v;
}

/* ---- Tier 3 ---- */

ML_API double ml_harmonic_number(long long n) {
    /* H_n = psi(n+1) + gamma (reuses the staged digamma core). */
    static const long double egamma = 0.57721566490153286060651209008240243L;
    long double h;
    double r;
    if (n < 1 || n > 1000000000000LL) return ml_make_nan();
    h = (long double)ml_digamma((double)n + 1.0) + egamma;
    r = (double)h;
    return ml_isfinite(r) ? r : ml_make_nan();
}

ML_API double ml_basel_partial(long long n) {
    /* Sum_{k<=n} 1/k^2, Kahan in LD; limit pi^2/6 (documented, not returned). */
    long double s = 0.0L, c = 0.0L;
    long long k;
    double r;
    if (n < 1 || n > 1000000LL) return ml_make_nan();
    for (k = 1; k <= n; k++) {
        long double t = 1.0L / ((long double)k * (long double)k);
        long double y = t - c;
        long double tt = s + y;
        c = (tt - s) - y;
        s = tt;
    }
    r = (double)s;
    return ml_isfinite(r) ? r : ml_make_nan();
}

ML_API double ml_faulhaber_sum(long long n, int p) {
    /* S_p via Bernoulli: 1/(p+1) sum C(p+1,j) B_j n^{p+1-j}.
     * LD evaluation (~1e-15 relative for n <= 1e6, p <= 10);
     * regression-guard tier, not exact-or-NaN (see header). */
    static const long double B[31] = {
        1.0L, -0.5L, 1.0L/6.0L, 0.0L, -1.0L/30.0L, 0.0L, 1.0L/42.0L, 0.0L,
        -1.0L/30.0L, 0.0L, 5.0L/66.0L, 0.0L, -691.0L/2730.0L, 0.0L,
        7.0L/6.0L, 0.0L, -3617.0L/510.0L, 0.0L, 43867.0L/798.0L, 0.0L,
        -174611.0L/330.0L, 0.0L, 854513.0L/138.0L, 0.0L,
        -236364091.0L/2730.0L, 0.0L, 8553103.0L/6.0L, 0.0L,
        -23749461029.0L/870.0L, 0.0L, 8615841276005.0L/14322.0L
    };
    long double nn;
    long double acc = 0.0L;
    int j;
    double r;
    if (n < 1 || n > 1000000LL || p < 0 || p > 10) return ml_make_nan();
    nn = (long double)n;
    for (j = 0; j <= p; j++) {
        /* C(p+1,j) exact in LD for p <= 10. (-1)^j pairs with the B1=-1/2
         * convention (all other odd B vanish); without it S_1 comes out
         * as (n^2-n)/2 instead of n(n+1)/2. */
        long double c = 1.0L;
        long double sgn = ((j & 1) != 0) ? -1.0L : 1.0L;
        int jj;
        int top = p + 1;
        int jj2 = j;
        int kk = top - j;
        int t;
        long double pw = 1.0L;
        int e;
        if (jj2 > kk) { t = jj2; jj2 = kk; kk = t; }
        for (jj = 0; jj < jj2; jj++) c = c * (long double)(top - jj) / (long double)(jj + 1);
        for (e = 0; e < p + 1 - j; e++) pw *= nn;
        acc += sgn * c * B[j] * pw;
    }
    acc /= (long double)(p + 1);
    r = (double)acc;
    return ml_isfinite(r) ? r : ml_make_nan();
}

ML_API ml_status_t ml_wythoff_cold(uint64_t k, uint64_t *a, uint64_t *b) {
    /* k-th cold: (floor(k phi), +k), phi = (1+sqrt5)/2. Exact for k<=9e15. */
    static const long double phi = 1.61803398874989484820458683436563811L;
    long double fa;
    uint64_t aa;
    if (ML_UNLIKELY(!a || !b || k < 1 || k > 9007199254740992ULL)) return ML_ERR_INVALID_ARG;
    fa = (long double)k * phi;
    if (!(fa >= 0.0L) || fa > 18446744073709551615.0L) return ML_ERR_OVERFLOW;
    aa = (uint64_t)fa;
    if (aa > UINT64_MAX - k) return ML_ERR_OVERFLOW;
    *a = aa;
    *b = aa + k;
    return ML_SUCCESS;
}

ML_API int ml_subtract_game_grundy(uint64_t n, const int *moves, int nmoves) {
    /* Subtraction-game Grundy via DP to n (cap 4096 => 16KB frame). */
    int g[4097];
    uint64_t i;
    int j;
    if (n > 4096 || nmoves < 1 || nmoves > 64 || !moves) return -1;
    for (j = 0; j < nmoves; j++) {
        if (moves[j] < 1 || moves[j] > 4096) return -1;
    }
    for (i = 0; i <= 4096; i++) g[i] = 0;
    for (i = 1; i <= n; i++) {
        int seen[65];
        int mex = 0;
        for (j = 0; j < 65; j++) seen[j] = 0;
        for (j = 0; j < nmoves; j++) {
            uint64_t mv = (uint64_t)moves[j];
            if (mv <= i) {
                int gv = g[(size_t)(i - mv)];
                if (gv >= 0 && gv < 65) seen[gv] = 1;
            }
        }
        while (mex < 65 && seen[mex] != 0) mex++;
        g[i] = mex;
    }
    return g[n];
}

ML_API int ml_schur_deg1(double a, double b, double c) {
    /* Schur r=1: sum a^3 + 3abc >= sum_sym a^2 b, a,b,c >= 0. */
    long double x;
    long double y;
    long double z;
    long double lhs;
    long double rhs;
    long double tol;
    if (ml_isnan(a) || ml_isnan(b) || ml_isnan(c)) return 0;
    if (!ml_isfinite(a) || !ml_isfinite(b) || !ml_isfinite(c)) return 0;
    if (a < 0.0 || b < 0.0 || c < 0.0) return 0;
    x = (long double)a;
    y = (long double)b;
    z = (long double)c;
    lhs = x*x*x + y*y*y + z*z*z + 3.0L*x*y*z;
    rhs = x*x*y + x*x*z + y*y*x + y*y*z + z*z*x + z*z*y;
    tol = 1e-12L * (1.0L + lhs + rhs);
    return (lhs + tol >= rhs) ? 1 : 0;
}

ML_API int ml_amgm_check(const double *data, int n) {
    /* Verifies AM >= GM on the instance (LD; always true for data >= 0). */
    long double s = 0.0L;
    long double ls = 0.0L;
    int i;
    long double am;
    long double gm;
    if (!data || n <= 0 || n > 1000000) return 0;
    for (i = 0; i < n; i++) {
        if (ml_isnan(data[i]) || !ml_isfinite(data[i]) || data[i] < 0.0) return 0;
        s += (long double)data[i];
        if (data[i] == 0.0) {
            /* GM = 0 <= AM: true unless NaN already handled. */
            long double am2 = s / (long double)n;
            return (am2 >= 0.0L) ? 1 : 0;
        }
        ls += __builtin_logl((long double)data[i]);
    }
    am = s / (long double)n;
    gm = __builtin_expl(ls / (long double)n);
    return (am + 1e-18L * (1.0L + am) >= gm) ? 1 : 0;
}

ML_API ml_status_t ml_vieta_cubic_sums(double a, double b, double c, double d,
                                       double *s1, double *s2, double *s3) {
    /* ax^3+bx^2+cx+d: s1=-b/a, s2=c/a, s3=-d/a. */
    if (ML_UNLIKELY(!s1 || !s2 || !s3)) return ML_ERR_INVALID_ARG;
    if (ML_UNLIKELY(ml_isnan(a) || ml_isnan(b) || ml_isnan(c) || ml_isnan(d))) return ML_ERR_NAN_INPUT;
    if (ML_UNLIKELY(!ml_isfinite(a) || !ml_isfinite(b) || !ml_isfinite(c) || !ml_isfinite(d))) return ML_ERR_NAN_INPUT;
    if (ML_UNLIKELY(a == 0.0)) return ML_ERR_INVALID_ARG;
    *s1 = -b / a;
    *s2 = c / a;
    *s3 = -d / a;
    if (ML_UNLIKELY(!ml_isfinite(*s1) || !ml_isfinite(*s2) || !ml_isfinite(*s3))) return ML_ERR_OVERFLOW;
    return ML_SUCCESS;
}

/* ---- Tier 0: recreational (AMC 8/10 kit) ---- */

ML_API uint64_t ml_digit_sum_u64(uint64_t n) {
    uint64_t s = 0;
    if (n == 0) return 0;
    while (n != 0) {
        s += n % 10ULL;
        n /= 10ULL;
    }
    return s;
}

ML_API uint64_t ml_digital_root(uint64_t n) {
    /* dr(n) = 1 + (n-1) mod 9 for n > 0; dr(0) = 0. */
    if (n == 0) return 0;
    return 1ULL + (n - 1ULL) % 9ULL;
}

ML_API int ml_is_palindrome_u64(uint64_t n) {
    uint64_t rev = 0;
    uint64_t t = n;
    while (t != 0) {
        uint64_t d = t % 10ULL;
        if (rev > (UINT64_MAX - d) / 10ULL) return 0;
        rev = rev * 10ULL + d;
        t /= 10ULL;
    }
    return (rev == n) ? 1 : 0;
}

ML_API uint64_t ml_collatz_steps(uint64_t n) {
    /* Total stopping time to 1. UINT64_MAX when the trajectory overflows
     * u64 or exceeds 10^6 steps (no claim on the open conjecture). */
    uint64_t steps = 0;
    if (n < 1) return UINT64_MAX;
    while (n != 1) {
        if (steps >= 1000000ULL) return UINT64_MAX;
        if ((n & 1ULL) == 0) {
            n /= 2ULL;
        } else {
            if (n > (UINT64_MAX - 1ULL) / 3ULL) return UINT64_MAX;
            n = 3ULL * n + 1ULL;
        }
        steps++;
    }
    return steps;
}

ML_API int ml_magic_square_check(const long long *a, int n) {
    /* Normal magic square: permutation of 1..n^2 with common line sums. */
    long long target;
    long long d1 = 0;
    long long d2 = 0;
    int seen[1025];
    int i;
    int j;
    if (!a || n < 3 || n > 32) return 0;
    for (i = 0; i <= n * n; i++) seen[i] = 0;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            long long v = a[(size_t)i * (size_t)n + (size_t)j];
            if (v < 1 || v > (long long)n * (long long)n) return 0;
            if (seen[v]) return 0;
            seen[v] = 1;
        }
    }
    target = 0;
    for (j = 0; j < n; j++) target += a[j];
    for (i = 0; i < n; i++) {
        long long rs = 0;
        long long cs = 0;
        for (j = 0; j < n; j++) {
            rs += a[(size_t)i * (size_t)n + (size_t)j];
            cs += a[(size_t)j * (size_t)n + (size_t)i];
        }
        if (rs != target || cs != target) return 0;
        d1 += a[(size_t)i * (size_t)n + (size_t)i];
        d2 += a[(size_t)i * (size_t)n + (size_t)(n - 1 - i)];
    }
    return (d1 == target && d2 == target) ? 1 : 0;
}
