#ifndef LIBMATHC_BITWISE_FP_H
#define LIBMATHC_BITWISE_FP_H
#include <stdint.h>
#include <string.h>
#include "ml_core.h"



#define ML_FP_SIGN_MASK  0x8000000000000000ULL
#define ML_FP_EXP_MASK   0x7FF0000000000000ULL
#define ML_FP_MANT_MASK  0x000FFFFFFFFFFFFFULL

// Pure bitwise classification without <math.h>
static inline int ml_fp_classify(double x) {
    uint64_t bits;
    memcpy(&bits, &x, sizeof(uint64_t)); // Safe extraction: prevents sNaN hardware traps & strict aliasing UB
    uint64_t exp = bits & ML_FP_EXP_MASK;
    uint64_t mant = bits & ML_FP_MANT_MASK;

    if (exp == ML_FP_EXP_MASK) return mant ? 4 : 3; // 4: NaN, 3: Inf
    if (exp == 0) return mant ? 1 : 0;              // 1: Subnormal, 0: Zero
    return 2;                                       // 2: Normal
}

static inline int ml_is_subnormal(double x) { return ml_fp_classify(x) == 1; }
static inline int ml_is_nan(double x) { return ml_fp_classify(x) == 4; }
static inline int ml_is_inf(double x) { return ml_fp_classify(x) == 3; }
#endif
