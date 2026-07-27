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
