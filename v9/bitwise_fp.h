#ifndef LIBMATHC_BITWISE_FP_H
#define LIBMATHC_BITWISE_FP_H
#include <stdint.h>

typedef union { double d; uint64_t u; } ml_fp_cast;

#define ML_FP_SIGN_MASK  0x8000000000000000ULL
#define ML_FP_EXP_MASK   0x7FF0000000000000ULL
#define ML_FP_MANT_MASK  0x000FFFFFFFFFFFFFULL

// Pure bitwise classification without <math.h>
inline int ml_fp_classify(double x) {
    ml_fp_cast c; c.d = x;
    uint64_t exp = c.u & ML_FP_EXP_MASK;
    uint64_t mant = c.u & ML_FP_MANT_MASK;

    if (exp == ML_FP_EXP_MASK) return mant ? 4 : 3; // 4: NaN, 3: Inf
    if (exp == 0) return mant ? 1 : 0;              // 1: Subnormal, 0: Zero
    return 2;                                       // 2: Normal
}

inline int ml_is_subnormal(double x) { return ml_fp_classify(x) == 1; }
inline int ml_is_nan(double x) { return ml_fp_classify(x) == 4; }
inline int ml_is_inf(double x) { return ml_fp_classify(x) == 3; }
#endif
