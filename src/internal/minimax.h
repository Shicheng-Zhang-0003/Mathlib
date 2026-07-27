#ifndef LIBMATHC_MINIMAX_H
#define LIBMATHC_MINIMAX_H

#include "ml_core.h"
#include "internal/payne_hanek.h"

// 19th-degree Maclaurin (Taylor) polynomial for sin(x) on [-pi/4, pi/4]
// Note: These are the exact Taylor series coefficients (1/k!). They achieve < 1 ULP
// on this interval due to the high degree. A true Remez minimax polynomial would
// have slightly different coefficients to minimize the maximum error uniformly.
static const double maclaurin_sin_coeffs[] = {
    1.0,
    -0.16666666666666666,
    0.008333333333333333,
    -0.0001984126984126984,
    2.7557319223985893e-06,
    -2.505210838544172e-08,
    1.6059043836821613e-10,
    -7.647163731819816e-13,
    2.811457254345521e-15,
    -8.220635816560923e-18
};

static inline double ml_minimax_sin_raw(double x) {
    double x2 = x * x;
    double result = maclaurin_sin_coeffs[9];
    for (int i = 8; i >= 0; i--) result = ML_FMA(result, x2, maclaurin_sin_coeffs[i]);
    return x * result;
}

// 18th-degree Maclaurin (Taylor) polynomial for cos(x) on [-pi/4, pi/4]
static const double maclaurin_cos_coeffs[] = {
    1.0,
    -0.5,
    0.041666666666666664,
    -0.001388888888888889,
    2.48015873015873e-05,
    -2.755731922398589e-07,
    2.08767569878681e-09,
    -1.1470745597729725e-11,
    4.779477332387385e-14,
    -1.5619206967218455e-16
};

static inline double ml_minimax_cos_raw(double x) {
    double x2 = x * x;
    double result = maclaurin_cos_coeffs[9];
    for (int i = 8; i >= 0; i--) result = ML_FMA(result, x2, maclaurin_cos_coeffs[i]);
    return result;
}

static inline double ml_minimax_sin(double x) {
    double y;
    int n = ml_rem_pio2(x, &y);
    if (ml_isnan(y)) return ml_make_nan();

    switch (n) {
        case 0: return  ml_minimax_sin_raw(y);
        case 1: return  ml_minimax_cos_raw(y);
        case 2: return -ml_minimax_sin_raw(y);
        case 3: return -ml_minimax_cos_raw(y);
    }
    return ml_make_nan();
}

static inline double ml_minimax_cos(double x) {
    double y;
    int n = ml_rem_pio2(x, &y);
    if (ml_isnan(y)) return ml_make_nan();

    switch (n) {
        case 0: return  ml_minimax_cos_raw(y);
        case 1: return -ml_minimax_sin_raw(y);
        case 2: return -ml_minimax_cos_raw(y);
        case 3: return  ml_minimax_sin_raw(y);
    }
    return ml_make_nan();
}

#endif
