#ifndef LIBMATHC_EXPONENTIAL_H
#define LIBMATHC_EXPONENTIAL_H

//Library header file for exponential and logarithmic functions
#include <math.h>
#ifndef M_E
#define M_E 2.71828182845904523536
#endif
#define math_e M_E
inline float exponential (float x);
inline float exponential (float x) {
    int k = 0;
    float y = x;
    while (y > 1) {y /= 2; k++;}
    while (y < -1) {y /= 2; k++;}
    float result = 1;
    float term = 1;
    for (int step = 1; step < 50; step++) {
        term *= y / step;
        result += term;
    }
    for (int i = 0; i < k; i++) {result *= result;}
    return result;
} inline float logarithm (float x);
inline float logarithm (float x) {
    if (x <= 0) {return 0.0 / 0.0;}
    int k = 0;
    float y = x;
    while (y > 1.5) {y /= 2; k++;}
    while (y < 0.5) {y *= 2; k--;}
    float z = y - 1;
    float result = z;
    for (int step = 2; step < 100; step++) {
        float additive = pow (z, (float) (step)) / step;
        if (step % 2 == 0) {result -= additive;}
        else {result += additive;}
    } return result + k * 0.69314718056;
} inline float power (float x, float y);
inline float power (float x, float y) {return exponential (y * logarithm (x));}
inline float logarithm_base (float x, float b);
inline float logarithm_base (float x, float b) {return logarithm (x) / logarithm (b);}
inline float hyperbolic_sine (float x);
inline float hyperbolic_sine (float x) {
    float result = x;
    for (int step = 3; step < 100; step += 2) {
        float additive = pow (x, (float) (step)) / (tgamma (step + 1));
        result += additive;
    } return result;
} inline float hyperbolic_cosine (float x);
inline float hyperbolic_cosine (float x) {
    float result = 1;
    for (int step = 2; step < 100; step += 2) {
        float additive = pow (x, (float) (step)) / (tgamma (step + 1));
        result += additive;
    } return result;
} inline float hyperbolic_tangent (float x);
inline float hyperbolic_tangent (float x) {return hyperbolic_sine (x) / hyperbolic_cosine (x);}
inline float inverse_hyperbolic_sine (float x);
inline float inverse_hyperbolic_sine (float x) {return logarithm (x + sqrt (x * x + 1));}
inline float inverse_hyperbolic_cosine (float x);
inline float inverse_hyperbolic_cosine (float x) {
    if (x < 1) {return 0.0 / 0.0;}
    return logarithm (x + sqrt (x * x - 1));
} inline float inverse_hyperbolic_tangent (float x);
inline float inverse_hyperbolic_tangent (float x) {
    if (x <= -1 || x >= 1) {return 0.0 / 0.0;}
    return 0.5 * logarithm ((1 + x) / (1 - x));
}

#endif
