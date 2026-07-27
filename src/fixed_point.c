#include "ml_compiler.h"
#include "ml_fixed_point.h"
#include <stdint.h>
#include <stddef.h>

/* v11S CLOSURE IP-17: fixed-point CORDIC portability hardening */

/*
 * CORDIC arctangent table for Q16.16.
 *
 * These constants are retained unchanged from the original v11S module
 * to preserve numerical behavior and existing test expectations.
 */
static const ml_q16_16_t fixed_cordic_atan[] = {
    51471, 30385, 16061, 8152, 4093, 2048, 1024, 512,
    256, 128, 64, 32, 16, 8, 4, 2, 1
};

/*
 * Defined arithmetic right shift for signed 64-bit integers.
 *
 * C99 leaves right-shift of negative signed integers implementation-defined.
 * This helper emulates the usual arithmetic-shift behavior portably:
 * floor(value / 2^shift).
 */
static inline int64_t ml_i64_asr(int64_t value, int shift) {
    if (shift == 0) {
        return value;
    }

    if (value >= 0) {
        return value >> shift;
    }

    int64_t divisor = (int64_t)1 << shift;
    int64_t pos = -value;

    /*
     * For negative values:
     * arithmetic shift == floor division by 2^shift
     */
    int64_t rounded = (pos + divisor - 1) >> shift;
    return -rounded;
}

/*
 * Defined arithmetic right shift for signed 32-bit integers.
 *
 * This is used by the CORDIC iteration to replace raw `x >> i`
 * and `y >> i`, which were implementation-defined for negative values.
 */
static inline ml_q16_16_t ml_q16_asr(ml_q16_16_t value, int shift) {
    if (shift == 0) {
        return value;
    }

    int64_t v = (int64_t)value;

    if (v >= 0) {
        return (ml_q16_16_t)(v >> shift);
    }

    int64_t divisor = (int64_t)1 << shift;
    int64_t pos = -v;
    int64_t rounded = (pos + divisor - 1) >> shift;

    return (ml_q16_16_t)(-rounded);
}

ML_API ml_q16_16_t ml_fixed_mul(ml_q16_16_t a, ml_q16_16_t b) {
    int64_t product = (int64_t)a * (int64_t)b;

    /*
     * Portable arithmetic shift instead of raw `>> 16`.
     */
    int64_t res = ml_i64_asr(product, 16);

    if (res > (int64_t)INT32_MAX) {
        return (ml_q16_16_t)INT32_MAX;
    }

    if (res < (int64_t)INT32_MIN) {
        return (ml_q16_16_t)INT32_MIN;
    }

    return (ml_q16_16_t)res;
}

ML_API void ml_cordic_sincos_fixed(ml_q16_16_t theta, ml_q16_16_t *sin_out, ml_q16_16_t *cos_out) {
    if (ML_UNLIKELY(sin_out == NULL || cos_out == NULL)) {
        return;
    }

    /*
     * Range reduction into [-pi, pi].
     *
     * The modulo result is then explicitly wrapped into the CORDIC
     * convergence domain.
     */
    theta = theta % ML_FIXED_TWO_PI;

    if (theta > ML_FIXED_PI) {
        theta -= ML_FIXED_TWO_PI;
    }

    if (theta < -ML_FIXED_PI) {
        theta += ML_FIXED_TWO_PI;
    }

    /*
     * Quadrant mapping:
     * CORDIC converges only inside [-pi/2, pi/2].
     */
    int negate_cos = 0;

    if (theta > ML_FIXED_HALF_PI) {
        theta = ML_FIXED_PI - theta;
        negate_cos = 1;
    } else if (theta < -ML_FIXED_HALF_PI) {
        theta = -ML_FIXED_PI - theta;
        negate_cos = 1;
    }

    ml_q16_16_t x = ML_FIXED_CORDIC_GAIN;
    ml_q16_16_t y = 0;
    ml_q16_16_t z = theta;

    /*
     * 16-iteration CORDIC kernel.
     *
     * The gain constant matches the 16-iteration regime.
     */
    for (int i = 0; i < 16; i++) {
        ml_q16_16_t x_shift = ml_q16_asr(x, i);
        ml_q16_16_t y_shift = ml_q16_asr(y, i);

        ml_q16_16_t x_new;
        ml_q16_16_t y_new;

        if (z >= 0) {
            x_new = x - y_shift;
            y_new = y + x_shift;
            z -= fixed_cordic_atan[i];
        } else {
            x_new = x + y_shift;
            y_new = y - x_shift;
            z += fixed_cordic_atan[i];
        }

        x = x_new;
        y = y_new;
    }

    if (negate_cos) {
        x = -x;
    }

    *cos_out = x;
    *sin_out = y;
}
