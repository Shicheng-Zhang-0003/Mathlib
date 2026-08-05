#ifndef LIBMATHC_PAYNE_HANEK_H
#define LIBMATHC_PAYNE_HANEK_H
/* MATHLIB_V12A1_PAYNE_HANEK_V6 */
/*
* RANGE REDUCTION FOR TRIGONOMETRIC FUNCTIONS
*
* Path 1  |x| <= 1e6 : 2-term Cody-Waite with error-free transforms.
* Path 2  |x| >  1e6 : Payne-Hanek with per-term integer/fractional
*                       separation and full-table sweep.
*
* V5 failed because it accumulated the full product x*(2/pi) in
* double precision without reducing mod 8, so the integer part
* (billions) swallowed the fractional part.
*
* V6 restores the V3 per-term approach:
*   - For each table term, extract integer contribution mod 4
*   - Kahan-accumulate only the fractional parts
*   - Full-table sweep (k=0..65) so no bits are dropped
*/
#include <string.h>
#include "ml_core.h"
#include "internal/error_free.h"

/* ---- Cody-Waite constants (path 1) ---- */
static const double
    ML_PH_PIO2_HI = 0x1.921fb54442d18p+0,
    ML_PH_PIO2_LO = 0x1.1a62633145c07p-54;

/* ---- 2/pi table, 24-bit chunks (Cephes / musl) ---- */
static const int32_t ml_two_over_pi[66] = {
    0xA2F983, 0x6E4E44, 0x1529FC, 0x2757D1, 0xF534DD, 0xC0DB62,
    0x95993C, 0x439041, 0xFE5163, 0xABDEBB, 0xC561B7, 0x246E3A,
    0x424DD2, 0xE00649, 0x2EEA09, 0xD1921C, 0xFE1DEB, 0x1CB129,
    0xA73EE8, 0x8235F5, 0x2EBB44, 0x84E99C, 0x7026B4, 0x5F7E41,
    0x3991D6, 0x398353, 0x39F49C, 0x845F8B, 0xBDF928, 0x3B1FF8,
    0x97FFDE, 0x05980F, 0xEF2F11, 0x8B5A0A, 0x6D1F6D, 0x367ECF,
    0x27CB09, 0xB74F46, 0x3F669E, 0x5FEA2D, 0x7527BA, 0xC7EBE5,
    0xF17B3D, 0x0739F7, 0x8A5292, 0xEA6BFB, 0x5FB11F, 0x8D5D08,
    0x560330, 0x46FC7B, 0x6BABF0, 0xCFBC20, 0x9AF436, 0x1DA9E3,
    0x91615E, 0xE61B08, 0x659985, 0x5F14A0, 0x68408D, 0xFFD880,
    0x4D7327, 0x310606, 0x1556CA, 0x73A8C9, 0x60E27B, 0xC08C6B
};

/* pi/2 high / low for reconstruction */
static const double
    ML_PH_PI2_HI = 0x1.921fb54442d18p+0,
    ML_PH_PI2_LO = 0x1.1a62633145c07p-54;
static const double ML_PH_TWO_OVER_PI = 0.63661977236758134308;

/*
* Process one term: prod * 2^shift.
*
* prod is an exact integer < 2^53.
* Extracts integer contribution mod 4 into *n.
* Kahan-accumulates fractional part into *acc / *comp.
*/
static inline void ml_ph_process_term(
    double prod, int shift, int *n, double *acc, double *comp
) {
    if (shift >= 2) {
        /* prod * 2^shift is integer divisible by 4. No contribution. */
        return;
    }
    if (shift == 1) {
        /* (prod*2) mod 4 = 2*(prod mod 2) */
        double fl = (double)((long long)prod & 1LL);
        *n = (*n + (int)(2.0 * fl)) & 3;
        return;
    }
    if (shift == 0) {
        /* prod mod 4 */
        double fl = (double)((long long)prod & 3LL);
        *n = (*n + (int)fl) & 3;
        return;
    }
    /* shift < 0: fractional term */
    double t = ml_ldexp_pure(prod, shift);
    double int_part = 0.0;
    double frac_part;
    if (shift >= -52) {
        /* t may have an integer part */
        int_part = (double)(long long)t;
        frac_part = t - int_part;
        *n = (*n + ((int)(long long)int_part & 3)) & 3;
    } else {
        /* t < 1, purely fractional */
        frac_part = t;
    }
    /* Kahan summation of fractional part */
    double y_k = frac_part - *comp;
    double t_k = *acc + y_k;
    *comp = (t_k - *acc) - y_k;
    *acc = t_k;
}

/*
* Payne-Hanek reduction for large arguments (|x| > 1e6).
*
* Computes:
*   *y = x mod (pi/2), in [-pi/4, pi/4]
*   returns quadrant n in {0, 1, 2, 3}
*/
static inline int ml_rem_pio2_large(double x, double *y) {
    uint64_t bits;
    double   ax   = ml_fabs(x);
    int      sign = (x < 0.0);

    memcpy(&bits, &ax, sizeof(uint64_t));
    int biased_e = (int)((bits >> 52) & 0x7FF);
    /* E: exponent such that ax = m * 2^E (m is 53-bit integer) */
    int E = biased_e - 1075;
    uint64_t m = (bits & 0x000FFFFFFFFFFFFFULL) | (1ULL << 52);

    /*
    * Split m into high 28 bits and low 25 bits.
    * m_hi has bits 52..25 (28 bits), m_lo has bits 24..0 (25 bits).
    */
    double m_hi = (double)(m >> 25);
    double m_lo = (double)(m & 0x1FFFFFFULL);

    /*
    * Full-table sweep.
    *
    * The skip logic inside ml_ph_process_term() handles terms that
    * don't contribute (shift >= 2). Processing all 66 entries
    * ensures no fractional bits are dropped for any |x|.
    */
    int k_start = 0;
    int k_end   = 65;

    /* Accumulate quadrant and fractional part */
    int    n         = 0;
    double frac_acc  = 0.0;
    double frac_comp = 0.0;

    for (int k = k_start; k <= k_end; k++) {
        double tk = (double)ml_two_over_pi[k];
        int base_shift = E - 24 * k - 24;

        /* m_hi term: m_hi * tk * 2^(base_shift + 25) */
        double prod_hi = m_hi * tk;  /* exact: 28+24 = 52 bits */
        ml_ph_process_term(prod_hi, base_shift + 25,
                           &n, &frac_acc, &frac_comp);

        /* m_lo term: m_lo * tk * 2^base_shift */
        double prod_lo = m_lo * tk;  /* exact: 25+24 = 49 bits */
        ml_ph_process_term(prod_lo, base_shift,
                           &n, &frac_acc, &frac_comp);
    }

    /* Extract integer part from fractional accumulator */
    int extra = (int)(long long)frac_acc;
    n = (n + (extra & 3)) & 3;
    double frac = frac_acc - (double)extra;

    /* Center fractional part in [-0.5, 0.5] */
    if (frac >  0.5) { frac -= 1.0; n = (n + 1) & 3; }
    if (frac < -0.5) { frac += 1.0; n = (n + 3) & 3; }

    /* Reconstruct: reduced_arg = frac * pi/2 */
    double result = ML_FMA(frac, ML_PH_PI2_HI, frac * ML_PH_PI2_LO);

    if (sign) {
        result = -result;
        n = (4 - n) & 3;
    }
    *y = result;
    return n;
}

/* ==================================================================
* UNIFIED ENTRY POINT
* ================================================================== */
static inline int ml_rem_pio2(double x, double *y) {
    if (ml_isnan(x) || ml_isinf(x)) {
        *y = ml_make_nan();
        return 0;
    }
    double ax = ml_fabs(x);

    /*
    * Small/medium arguments: fast Cody-Waite.
    * Accurate to < 1 ULP for |x| <= 1e6.
    */
    if (ax <= 1.0e6) {
        double fn = ml_round(x * ML_PH_TWO_OVER_PI);
        long long n_ll = (long long)fn;
        int n = (int)(n_ll % 4);
        if (n < 0) n += 4;

        double p     = fn * ML_PH_PIO2_HI;
        double p_err = ML_FMA(fn, ML_PH_PIO2_HI, -p);
        double r1, r1_err;
        r1 = ml_two_sum(x, -p, &r1_err);
        double r2 = r1_err - p_err - (fn * ML_PH_PIO2_LO);
        *y = r1 + r2;
        return n;
    }

    /*
    * Large arguments: Payne-Hanek.
    * Works for the full double range up to ~1.8e308.
    */
    return ml_rem_pio2_large(x, y);
}

#endif /* LIBMATHC_PAYNE_HANEK_H */
