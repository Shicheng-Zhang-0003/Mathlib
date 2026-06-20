#ifndef LIBMATHC_INTEGRAL_H
#define LIBMATHC_INTEGRAL_H

#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_E
#define M_E 2.71828182845904523536
#endif
#define math_pi M_PI
#define math_e M_E
inline float factorial_float (float x);
inline float factorial_float (float x) {return (sqrt (2 * math_pi * x)) * (pow ((x / math_e), x));}
inline float integral_traditional (float a, float b, float exponent, float additive, float d) {
    float result = 0.0;
    float x = a;
    while (x < b) {result += (pow (x, exponent) + additive) * d;
    x += d;} return result;
} inline float gamma_new (float x);
inline float gamma_new (float x) {
    if (x <= 0) {return 0.0 / 0.0;}
    if (x > 2) {return (x - 1) * gamma_new (x - 1);}
    if (x < 1) {return gamma_new (x + 1) / x;}
    float z = x - 1;
    float p = -0.193527818 + z * 0.035868343;
    p = 0.482199394 + z * p;
    p = -0.756704078 + z * p;
    p = 0.918206857 + z * p;
    p = -0.897056937 + z * p;
    p = 0.989028236 + z * p;
    p = -0.577191652 + z * p;
    return 1 + z * p;
}

#endif
