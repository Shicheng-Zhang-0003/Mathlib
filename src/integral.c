#include "ml_compiler.h"
#include "ml_integral.h"

/* v11S CLOSURE IP-14: integral / gamma semantics cleanup */

/*
 * v11S contract:
 * ml_factorial_float(x) means x! = Gamma(x + 1).
 *
 * This is the semantically correct interpretation users expect.
 * Negative finite inputs are invalid and return NaN.
 */
ML_API double ml_factorial_float(double x) {
    if (ml_isnan(x)) return x;
    if (x < 0.0) return ml_make_nan();
    if (ml_isinf(x)) return ml_make_inf(0);
    if (x == 0.0) return 1.0;

    return ml_gamma_new(x + 1.0);
}

/*
 * Educational / experimental numerical integrator.
 *
 * This is not a primary high-accuracy quadrature API.
 * It is retained for compatibility and simple classroom-style integration.
 *
 * v11S hardening:
 * - rejects non-finite parameters
 * - keeps O(1) temporal ceiling
 * - detects floating-point stall (x + d == x)
 * - returns NaN instead of silently producing garbage
 */
ML_API double ml_integral_traditional(double a, double b, double exponent, double additive, double d) {
    if (ml_isnan(a) || ml_isnan(b) || ml_isnan(exponent) ||
        ml_isnan(additive) || ml_isnan(d)) {
        return ml_make_nan();
    }

    if (ml_isinf(a) || ml_isinf(b) || ml_isinf(exponent) ||
        ml_isinf(additive) || ml_isinf(d)) {
        return ml_make_nan();
    }

    if (d == 0.0) {
        return ml_make_nan();
    }

    if ((d > 0.0 && a >= b) || (d < 0.0 && a <= b)) {
        return 0.0;
    }

    double result = 0.0;
    double x = a;

    const int max_steps = 10000000; /* O(1) temporal ceiling */

    for (int step = 0; step < max_steps; step++) {
        if ((d > 0.0 && x >= b) || (d < 0.0 && x <= b)) {
            return result;
        }

        double term = ml_pow(x, exponent) + additive;

        if (ml_isnan(term)) {
            return ml_make_nan();
        }

        result += term * d;

        /*
         * Floating-point stall detection:
         * If d is so small relative to x that x + d == x, then the loop
         * can never terminate naturally. Abort with NaN.
         */
        double next_x = x + d;

        if (ML_UNLIKELY(next_x == x)) {
            return ml_make_nan();
        }

        x = next_x;
    }

    /* Step limit exceeded: signal failure instead of returning a wrong answer. */
    return ml_make_nan();
}

/*
 * v11S contract:
 * ml_gamma_new() currently supports positive real inputs only.
 *
 * Negative non-integer extension is deferred to v12.
 * For v11S closure, the important fixes are:
 * - sane NaN / Inf behavior
 * - exact small-integer behavior
 * - no weird recursion / stack behavior
 * - predictable overflow to +Inf
 */
ML_API double ml_gamma_new(double x) {
    if (ml_isnan(x)) return x;

    if (ml_isinf(x)) {
        return x > 0.0 ? ml_make_inf(0) : ml_make_nan();
    }

    if (x <= 0.0) {
        return ml_make_nan(); /* positive-domain-only contract */
    }

    if (x > 171.0) {
        return ml_make_inf(0); /* Gamma(172) overflows double */
    }

    /*
     * Iterative reduction:
     * Gamma(x + 1) = x * Gamma(x)
     *
     * This prevents deep recursion on large positive inputs.
     */
    double result = 1.0;

    while (x > 2.0) {
        x -= 1.0;
        result *= x;
    }

    /*
     * Exact anchors:
     * Gamma(1) = 1
     * Gamma(2) = 1
     *
     * This also makes small integer factorials exact after reduction.
     */
    if (x == 1.0 || x == 2.0) {
        return result;
    }

    if (x < 1.0) {
        /*
         * Bounded recursion:
         * x + 1 is in [1.0, 2.0), so the next call hits the base case
         * immediately. Maximum stack depth is strictly 2.
         */
        return result * ml_gamma_new(x + 1.0) / x;
    }

    /*
     * Base case: x is in (1.0, 2.0)
     *
     * Approximation polynomial for the [1,2] domain.
     * This is not a full production-grade lgamma / tgamma implementation,
     * but it is stable and adequate for v11S closure.
     */
    double z = x - 1.0;

    double p = -0.193527818 + z * 0.035868343;
    p = 0.482199394 + z * p;
    p = -0.756704078 + z * p;
    p = 0.918206857 + z * p;
    p = -0.897056937 + z * p;
    p = 0.989028236 + z * p;
    p = -0.577191652 + z * p;

    return result * (1.0 + z * p);
}
