#ifndef MATHLIB_PROFILES_H
#define MATHLIB_PROFILES_H

#include "ml_core.h"
#include "fast_math.h"

// v11A1: profiles.h no longer overrides core math functions via macros.
// Core math (ml_sin, ml_exp) is handled cleanly in ml_trig.h and ml_exp_log.h.
// This file now only configures explicit fast-math approximations.

#if defined(MATHLIB_PROFILE_GRAPHICS)
    // Graphics profile exposes fast approximations via explicit ml_fast_ calls
    #define ml_rsqrt(x) ml_fast_rsqrt(x)
#elif defined(MATHLIB_PROFILE_EMBEDDED)
    // Embedded hides FPU-dependent fast math entirely
    // (ml_fast_rsqrt is not defined, forcing a compile error if used)
#else
    // Scientific strictness
    #define ml_rsqrt(x) (1.0 / ml_sqrt(x))
#endif

#endif
