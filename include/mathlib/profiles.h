#ifndef MATHLIB_PROFILES_H
#define MATHLIB_PROFILES_H

#include "ml_compiler.h"
#include "fast_math.h"

/* FIX: Eradicate macro poisoning. Use inline functions to prevent double-evaluation bugs. */
#if defined(MATHLIB_PROFILE_GRAPHICS)
ML_INLINE double ml_rsqrt(double x) {
    return ml_fast_rsqrt(x);
}
#elif defined(MATHLIB_PROFILE_EMBEDDED)
/* Embedded hides FPU-dependent fast math entirely */
#else
ML_INLINE double ml_rsqrt(double x) {
    return 1.0 / ml_sqrt(x);
}
#endif

#endif /* MATHLIB_PROFILES_H */
