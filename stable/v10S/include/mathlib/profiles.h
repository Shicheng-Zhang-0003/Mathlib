#ifndef MATHLIB_PROFILES_H
#define MATHLIB_PROFILES_H

#include "ml_core.h"
#include "minimax.h"
#include "fixed_point.h"
#include "cordic.h"
#include "trigonometry.h"
#include "fast_math.h"

#if defined(MATHLIB_PROFILE_GRAPHICS)
    // Speed Demon: Minimax + Fast RSqrt
    #define ml_sin(x) ml_minimax_sin(x)
    #define ml_cos(x) ml_minimax_cos(x)
    #define ml_sqrt(x) (1.0 / ml_fast_rsqrt(x))
    #define ml_rsqrt(x) ml_fast_rsqrt(x)
#elif defined(MATHLIB_PROFILE_EMBEDDED)
    // Survivalist: Pure Shift-and-Add CORDIC
    static inline double __ml_cordic_sin(double x) { double s, c; ml_cordic_sincos(x, &s, &c); return s; }
    static inline double __ml_cordic_cos(double x) { double s, c; ml_cordic_sincos(x, &s, &c); return c; }
    #define ml_sin(x) __ml_cordic_sin(x)
    #define ml_cos(x) __ml_cordic_cos(x)
    #define ml_sqrt(x) (1.0 / ml_fast_rsqrt(x))
    #define ml_rsqrt(x) ml_fast_rsqrt(x)
#else
    // Scientific (Default): Strict IEEE-754
    static inline double __ml_sci_sin(double x) { return ml_minimax_sin(x); }
    static inline double __ml_sci_cos(double x) { return ml_minimax_cos(x); }
    static inline double __ml_sci_sqrt(double x) { return ml_sqrt(x); }
    static inline double __ml_sci_rsqrt(double x) { return 1.0 / ml_sqrt(x); }
    #define ml_sin(x) __ml_sci_sin(x)
    #define ml_cos(x) __ml_sci_cos(x)
    #define ml_sqrt(x) __ml_sci_sqrt(x)
    #define ml_rsqrt(x) __ml_sci_rsqrt(x)
#endif

#endif
