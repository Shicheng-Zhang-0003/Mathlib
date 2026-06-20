#ifndef LIBMATHC_POLYNOMIAL_H
#define LIBMATHC_POLYNOMIAL_H

//Library header file for polynomial operations
#include <math.h>
inline float polynomial_eval (float *coeffs, int degree, float x);
inline float polynomial_eval (float *coeffs, int degree, float x) {
    float result = coeffs[degree];
    for (int i = degree - 1; i >= 0; i--) {result = result * x + coeffs[i];}
    return result;
} inline void polynomial_derivative (float *coeffs, int degree, float *out);
inline void polynomial_derivative (float *coeffs, int degree, float *out) {
    for (int i = 0; i < degree; i++) {out[i] = coeffs[i + 1] * (i + 1);}
} inline float polynomial_newton (float *coeffs, int degree, float x0, float epsilon, int max_iter);
inline float polynomial_newton (float *coeffs, int degree, float x0, float epsilon, int max_iter) {
    float x = x0;
    for (int iter = 0; iter < max_iter; iter++) {
        float fx = coeffs[degree];
        float dfx = 0;
        for (int i = degree - 1; i >= 0; i--) {dfx = dfx * x + fx;
        fx = fx * x + coeffs[i];}
        if (fabsf (dfx) < epsilon) {return 0.0 / 0.0;}
        float x_next = x - fx / dfx;
        if (fabsf (x_next - x) < epsilon) {return x_next;}
        x = x_next;
    } return x;
}



#endif
