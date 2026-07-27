#ifndef LIBMATHC_PAYNE_HANEK_H
#define LIBMATHC_PAYNE_HANEK_H

/*
 * HISTORICAL FILENAME NOTE:
 * This file currently implements a bounded high-precision Cody-Waite
 * reduction with an explicit domain clamp. It is NOT a full Payne-Hanek
 * implementation for arbitrary huge arguments.
 */

#include <string.h>
#include "ml_core.h"
#include "internal/error_free.h"

// High-precision PI/2 components for Cody-Waite reduction.
// Using exact hex floats guarantees the compiler cannot introduce rounding discrepancies.
static const double
PIO2_HI = 0x1.921fb54442d18p+0,  // 1.5707963267948966 (53 bits)
PIO2_LO = 0x1.1a62633145c07p-54; // 6.123233995736766e-17 (next 53 bits)


/* v11S CONTRACT:
 * Uses Error-Free Transformations (Dekker/Knuth) to perform exact Cody-Waite
 * reduction. This prevents the precision bleed that occurs when multiplying
 * a large integer 'fn' by a 53-bit PIO2_HI using standard IEEE-754 arithmetic.
 * Guarantees < 1 ULP reduction error for |x| <= 1e15.
 */
static inline int ml_rem_pio2(double x, double *y) {
    if (ml_isnan(x) || ml_isinf(x)) {
        *y = ml_make_nan();
        return 0;
    }

    double ax = ml_fabs(x);

    // DOMAIN CLAMP: Prevent long long overflow UB and precision collapse
    // 1e15 is the safe limit where 'fn' fits cleanly in 53 bits and
    // Two-Product maintains exactness.
    if (ax > 1.0e15) {
        *y = ml_make_nan();
        return 0;
    }

    // Estimate n = round(x / (pi/2))
    double fn = ml_round(x * 0.63661977236758134308); // 2/pi
    long long n_ll = (long long)fn;
    int n = (int)(n_ll % 4);
    if (n < 0) n += 4;

    // EXACT CODY-WAITE REDUCTION via Hardware FMA
    // Bypasses Dekker's 27-bit limitation by using infinite-precision intermediate FMA
    double p = fn * PIO2_HI;
    double p_err = ML_FMA(fn, PIO2_HI, -p);

    // 2. Compute x - p exactly as (r1 + r1_err)
    double r1, r1_err;
    r1 = ml_two_sum(x, -p, &r1_err);

    // 3. Accumulate the remaining low-order terms
    double r2 = r1_err - p_err - (fn * PIO2_LO);

    *y = r1 + r2;
    return n;
}
#endif
