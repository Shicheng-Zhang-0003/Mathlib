#ifndef ML_IEEE_EXACT_H
#define ML_IEEE_EXACT_H

/* ============================================================================
 * MATHLIB v11S: EXACT IEEE-754 DECOMPOSITION / RECOMPOSITION LAYER
 *
 * PURPOSE:
 * Provides a single, authoritative model for converting between IEEE-754
 * double-precision bit patterns and exact integer-significand representations.
 *
 * All downstream core primitives (ml_ldexp_pure, ml_frexp_pure, ml_fmod)
 * MUST be rebuilt on this layer to eliminate ad-hoc bit-manipulation bugs,
 * subnormal flushing, and silent precision loss.
 *
 * INVARIANT:
 * For any finite nonzero double x:
 *   x == sign * sig * 2^exp   (exactly, in infinite-precision arithmetic)
 *
 * where sig is a 53-bit integer (bit 52 set for normals, < 2^52 for subnormals)
 * and exp is the unbiased power-of-two exponent.
 * ========================================================================== */

#include <stdint.h>
#include <string.h>
#include "ml_core.h"

/* Classification of IEEE-754 double-precision values */
typedef enum {
    ML_FP_ZERO      = 0,
    ML_FP_SUBNORMAL = 1,
    ML_FP_NORMAL    = 2,
    ML_FP_INF       = 3,
    ML_FP_NAN       = 4
} ml_fp_kind_t;

/* Exact decomposition of a double into integer significand + exponent */
typedef struct {
    uint64_t       sig;   /* Exact integer significand (53-bit for normals) */
    int            exp;   /* Unbiased exponent: value = sign * sig * 2^exp */
    int            sign;  /* 0 = positive, 1 = negative */
    ml_fp_kind_t   kind;  /* Classification */
} ml_fp_parts_t;

/* ============================================================================
 * DECOMPOSE: double -> (sig, exp, sign, kind)
 *
 * Extracts the exact integer significand and unbiased exponent from any
 * IEEE-754 double. No precision is lost. No subnormals are flushed.
 *
 * NORMALS:   sig has bit 52 set (hidden bit restored), exp = biased - 1023 - 52
 * SUBNORMALS: sig < 2^52 (no hidden bit), exp = -1074 (fixed)
 * ZERO:      sig = 0, exp = 0
 * INF/NAN:   sig = raw mantissa bits, exp = 0, kind distinguishes them
 * ========================================================================== */
static inline ml_fp_parts_t ml_fp_decompose(double x) {
    uint64_t bits;
    memcpy(&bits, &x, sizeof(uint64_t));

    ml_fp_parts_t p;
    p.sign = (int)(bits >> 63);

    uint64_t exp_bits = (bits >> 52) & 0x7FFULL;
    uint64_t mant     = bits & 0x000FFFFFFFFFFFFFULL;

    /* Special values: exponent field all ones */
    if (exp_bits == 0x7FFULL) {
        p.kind = mant ? ML_FP_NAN : ML_FP_INF;
        p.sig  = mant;
        p.exp  = 0;
        return p;
    }

    /* Zero or subnormal: exponent field all zeros */
    if (exp_bits == 0) {
        if (mant == 0) {
            p.kind = ML_FP_ZERO;
            p.sig  = 0;
            p.exp  = 0;
            return p;
        }
        /* Subnormal: no hidden bit, exponent is fixed at -1074 */
        p.kind = ML_FP_SUBNORMAL;
        p.sig  = mant;
        p.exp  = -1074;
        return p;
    }

    /* Normal: restore hidden bit 52 */
    p.kind = ML_FP_NORMAL;
    p.sig  = mant | (1ULL << 52);
    p.exp  = (int)exp_bits - 1023 - 52;
    return p;
}

/* ============================================================================
 * COMPOSE: (sig, exp, sign) -> double
 *
 * Reconstructs an IEEE-754 double from an exact integer significand and
 * unbiased exponent. Handles:
 *
 * - Normalization of unnormalized significands
 * - Correct subnormal generation (gradual underflow)
 * - Round-to-nearest-even on subnormal underflow
 * - Overflow to signed infinity
 * - Underflow to signed zero
 * - Signed zero preservation when sig == 0
 *
 * This function NEVER flushes subnormals to zero unless the true
 * mathematical result is smaller than the smallest representable subnormal.
 * ========================================================================== */
static inline double ml_fp_compose(uint64_t sig, int exp, int sign) {
    /* Zero significand -> signed zero */
    if (sig == 0) {
        uint64_t bits = sign ? 0x8000000000000000ULL : 0ULL;
        double d;
        memcpy(&d, &bits, sizeof(double));
        return d;
    }

    /* Step 1: Normalize so that bit 52 is set (if possible within range) */
    while (sig < (1ULL << 52) && exp > -1074) {
        sig <<= 1;
        exp--;
    }

    /* Step 2: Check if we have a normal-range result */
    if (sig >= (1ULL << 52)) {
        int biased = exp + 52 + 1023;  /* Convert to IEEE biased exponent */

        /* Overflow: exponent too large -> signed infinity */
        if (biased >= 0x7FF) {
            return ml_make_inf(sign);
        }

        /* Normal range: biased exponent in [1, 2046] */
        if (biased >= 1) {
            uint64_t bits =
                ((uint64_t)sign << 63) |
                ((uint64_t)biased << 52) |
                (sig & 0x000FFFFFFFFFFFFFULL);
            double d;
            memcpy(&d, &bits, sizeof(double));
            return d;
        }

        /* Step 3: Subnormal range (biased <= 0) */
        /* We need to shift the significand right to denormalize it */
        int shift = 1 - biased;  /* Number of bits to shift right */

        if (shift >= 64) {
            /* Complete underflow -> signed zero */
            uint64_t bits = sign ? 0x8000000000000000ULL : 0ULL;
            double d;
            memcpy(&d, &bits, sizeof(double));
            return d;
        }

        uint64_t mant = sig >> shift;

        /* Round-to-nearest-even on the dropped bits */
        if (shift > 0) {
            uint64_t dropped_mask = (1ULL << shift) - 1ULL;
            uint64_t dropped = sig & dropped_mask;
            uint64_t half = 1ULL << (shift - 1);

            if (dropped > half || (dropped == half && (mant & 1ULL))) {
                mant++;
                /* Rounding up may push us back into normal range */
                if (mant >= (1ULL << 52)) {
                    uint64_t bits =
                        ((uint64_t)sign << 63) |
                        (1ULL << 52);  /* biased exp = 1, mantissa = 0 */
                    double d;
                    memcpy(&d, &bits, sizeof(double));
                    return d;
                }
            }
        }

        /* Subnormal: biased exponent field = 0, mantissa = mant */
        uint64_t bits = ((uint64_t)sign << 63) | mant;
        double d;
        memcpy(&d, &bits, sizeof(double));
        return d;
    }

    /* Step 4: sig < 2^52 and exp == -1074 (subnormal that couldn't normalize) */
    if (exp < -1074) {
        /* Complete underflow -> signed zero */
        uint64_t bits = sign ? 0x8000000000000000ULL : 0ULL;
        double d;
        memcpy(&d, &bits, sizeof(double));
        return d;
    }

    /* Subnormal with sig already in position */
    uint64_t bits = ((uint64_t)sign << 63) | sig;
    double d;
    memcpy(&d, &bits, sizeof(double));
    return d;
}

#endif /* ML_IEEE_EXACT_H */
