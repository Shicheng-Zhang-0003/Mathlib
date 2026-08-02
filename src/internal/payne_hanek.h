#ifndef LIBMATHC_PAYNE_HANEK_H
#define LIBMATHC_PAYNE_HANEK_H

/* MATHLIB_V12A1_PAYNE_HANEK_V2 */
#include <string.h>
#include "ml_core.h"
#include "internal/error_free.h"

static const double
ML_PH_PIO2_HI = 0x1.921fb54442d18p+0,
ML_PH_PIO2_LO = 0x1.1a62633145c07p-54;

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

static const double
ML_PH_PI2_HI = 0x1.921fb54442d18p+0,
ML_PH_PI2_LO = 0x1.1a62633145c07p-54;
static const double ML_PH_TWO_OVER_PI = 0.63661977236758134308;

static inline int ml_rem_pio2_large(double x, double *y) {
    uint64_t bits;
    double ax = ml_fabs(x);
    int sign = (x < 0.0);
    memcpy(&bits, &ax, sizeof(uint64_t));
    int biased_e = (int)((bits >> 52) & 0x7FF);
    int E = biased_e - 1075;
    uint64_t m = (bits & 0x000FFFFFFFFFFFFFULL) | (1ULL << 52);
    double m_hi = (double)(m >> 25);
    double m_lo = (double)(m & 0x1FFFFFFULL);
    int k_start = (E - 77) / 24;
    if (k_start < 0) k_start = 0;
    int k_end = (E + 53) / 24;
    if (k_end > 65) k_end = 65;
    double q_hi = 0.0;
    double q_lo = 0.0;
    for (int k = k_start; k <= k_end; k++) {
        int shift = E - 24 * k - 24;
        double tk = (double)ml_two_over_pi[k];
        double prod_hi = m_hi * tk;
        double prod_lo = m_lo * tk;
        q_hi += ml_ldexp_pure(prod_hi, shift + 25);
        q_lo += ml_ldexp_pure(prod_lo, shift);
    }
    double q = q_hi + q_lo;
    double n_d = ml_round(q);
    long long n_ll = (long long)n_d;
    int n = (int)(n_ll % 4);
    if (n < 0) n += 4;
    double frac = q - n_d;
    if (frac > 0.5) { frac -= 1.0; n = (n + 1) & 3; }
    else if (frac < -0.5) { frac += 1.0; n = (n + 3) & 3; }
    double result = frac * ML_PH_PI2_HI + frac * ML_PH_PI2_LO;
    if (sign) { result = -result; n = (4 - n) & 3; }
    *y = result;
    return n;
}

static inline int ml_rem_pio2(double x, double *y) {
    if (ml_isnan(x) || ml_isinf(x)) { *y = ml_make_nan(); return 0; }
    double ax = ml_fabs(x);
    if (ax <= 1.0e6) {
        double fn = ml_round(x * ML_PH_TWO_OVER_PI);
        long long n_ll = (long long)fn;
        int n = (int)(n_ll % 4);
        if (n < 0) n += 4;
        double p = fn * ML_PH_PIO2_HI;
        double p_err = ML_FMA(fn, ML_PH_PIO2_HI, -p);
        double r1, r1_err;
        r1 = ml_two_sum(x, -p, &r1_err);
        double r2 = r1_err - p_err - (fn * ML_PH_PIO2_LO);
        *y = r1 + r2;
        return n;
    }
    return ml_rem_pio2_large(x, y);
}

#endif /* LIBMATHC_PAYNE_HANEK_H */
