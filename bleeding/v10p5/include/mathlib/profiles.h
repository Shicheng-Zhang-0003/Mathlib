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
        #define ml_rsqrt(x) ml_fast_rsqrt(x)
#elif defined(MATHLIB_PROFILE_EMBEDDED)
    // Survivalist: True Q16.16 Bare-Metal Fixed-Point CORDIC (Zero FPU, Zero Division)
    static inline double __ml_cordic_sin_fixed(double x) {
        ml_q16_16_t f_in = (ml_q16_16_t)(x * 65536.0);
        ml_q16_16_t s, c;
        ml_cordic_sincos_fixed(f_in, &s, &c);
        return (double)s / 65536.0;
    }
    static inline double __ml_cordic_cos_fixed(double x) {
        ml_q16_16_t f_in = (ml_q16_16_t)(x * 65536.0);
        ml_q16_16_t s, c;
        ml_cordic_sincos_fixed(f_in, &s, &c);
        return (double)c / 65536.0;
    }
    #define ml_sin(x) __ml_cordic_sin_fixed(x)
    #define ml_cos(x) __ml_cordic_cos_fixed(x)
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
