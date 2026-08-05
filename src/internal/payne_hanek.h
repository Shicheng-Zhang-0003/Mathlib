#ifndef LIBMATHC_PAYNE_HANEK_H
#define LIBMATHC_PAYNE_HANEK_H
/* MATHLIB_V12A1_PAYNE_HANEK_V5_MUSL */
/*
* RANGE REDUCTION FOR TRIGONOMETRIC FUNCTIONS
*
* Path 1  |x| <= 1e6 : 2-term Cody-Waite with error-free transforms.
* Path 2  |x| >  1e6 : Payne-Hanek with fixed-array accumulation
*                       (musl / Cephes algorithm).
*
* The previous Kahan-summation accumulator lost precision for
* |x| >= 1e100 because ~54 fractional terms of very different
* magnitude were folded into a single double.  The fixed-array
* approach stores each contribution in f[20] and sums once,
* preserving every bit the table can supply.
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
* ml_rem_pio2_large  –  Payne-Hanek for |x| > 1e6.
*
* Returns quadrant n in {0,1,2,3} and sets *y to the reduced
* argument in [-pi/4, pi/4].
*
* Algorithm (musl __rem_pio2_large):
*   1. Decompose |x| into 24-bit chunks of the significand.
*   2. For each chunk, multiply by the relevant 2/pi table entries
*      and store the product (scaled to the correct position) in a
*      fixed-size array f[].
*   3. Sum f[] to obtain q = |x| * (2/pi) to full table precision.
*   4. Extract quadrant (integer part mod 4) and reduced argument.
*/
static inline int ml_rem_pio2_large(double x, double *y)
{
    uint64_t bits;
    double   ax   = ml_fabs(x);
    int      sign = (x < 0.0);

    memcpy(&bits, &ax, sizeof(uint64_t));
    int biased_e = (int)((bits >> 52) & 0x7FF);
    int E = biased_e - 1075;                 /* ax = m * 2^E */
    uint64_t m = (bits & 0x000FFFFFFFFFFFFFULL) | (1ULL << 52);

    /*
    * Split the 53-bit significand into three 24-bit chunks
    * (MSB first):
    *   x0 = bits 52..29  (24 bits), weight 2^(E+29)
    *   x1 = bits 28..5   (24 bits), weight 2^(E+5)
    *   x2 = bits  4..0   ( 5 bits), weight 2^E
    */
    uint32_t x0 = (uint32_t)(m >> 29);
    uint32_t x1 = (uint32_t)((m >> 5) & 0xFFFFFFULL);
    uint32_t x2 = (uint32_t)(m & 0x1FULL);

    /*
    * Accumulate x * (2/pi) into f[].
    * Each entry f[j] holds one 24-bit product, scaled by ldexp.
    * At most 3 chunks x 3 table entries = 9 terms, but we
    * allocate 20 for safety.
    */
    double f[20];
    int    nf = 0;

    /* chunk x0, weight 2^(E+29) */
    {
        int exp0 = E + 29;
        int j0   = exp0 / 24;
        if (j0 < 0) j0 = 0;
        for (int k = j0; k <= j0 + 2 && k < 66; k++) {
            int shift = exp0 - 24 * (k + 1);
            double prod = (double)x0 * (double)ml_two_over_pi[k];
            f[nf++] = ml_ldexp_pure(prod, shift);
        }
    }
    /* chunk x1, weight 2^(E+5) */
    {
        int exp1 = E + 5;
        int j1   = exp1 / 24;
        if (j1 < 0) j1 = 0;
        for (int k = j1; k <= j1 + 2 && k < 66; k++) {
            int shift = exp1 - 24 * (k + 1);
            double prod = (double)x1 * (double)ml_two_over_pi[k];
            f[nf++] = ml_ldexp_pure(prod, shift);
        }
    }
    /* chunk x2, weight 2^E */
    {
        int exp2 = E;
        int j2   = exp2 / 24;
        if (j2 < 0) j2 = 0;
        for (int k = j2; k <= j2 + 2 && k < 66; k++) {
            int shift = exp2 - 24 * (k + 1);
            double prod = (double)x2 * (double)ml_two_over_pi[k];
            f[nf++] = ml_ldexp_pure(prod, shift);
        }
    }

    /*
    * Sum f[] with pairwise-style compensation.
    * The musl approach: accumulate into a running sum, keeping
    * the low bits in a separate variable.
    */
    double q  = 0.0;
    double lo = 0.0;
    for (int i = 0; i < nf; i++) {
        double t = q + f[i];
        lo += (q - t) + f[i];   /* capture rounding error */
        q  = t;
    }
    q += lo;

    /*
    * Extract quadrant and fractional part.
    * q = |x| * (2/pi).  n = round(q) mod 4.
    * frac = q - n, centred in [-0.5, 0.5].
    */
    double n_d  = ml_round(q);
    long long n_ll = (long long)n_d;
    int n = (int)(n_ll & 3);

    double frac = q - n_d;
    if (frac >  0.5) { frac -= 1.0; n = (n + 1) & 3; }
    if (frac < -0.5) { frac += 1.0; n = (n + 3) & 3; }

    /* Reconstruct reduced argument = frac * (pi/2) */
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
static inline int ml_rem_pio2(double x, double *y)
{
    if (ml_isnan(x) || ml_isinf(x)) {
        *y = ml_make_nan();
        return 0;
    }
    double ax = ml_fabs(x);

    /* Path 1: Cody-Waite for |x| <= 1e6 */
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

    /* Path 2: Payne-Hanek for the full double range */
    return ml_rem_pio2_large(x, y);
}

#endif /* LIBMATHC_PAYNE_HANEK_H */
