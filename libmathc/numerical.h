#ifndef LIBMATHC_NUMERICAL_H
#define LIBMATHC_NUMERICAL_H

//Library header file for numerical methods
#include <math.h>
inline float newton_raphson (float (*f)(float), float (*df)(float), float x0, float epsilon, int max_iter);
inline float newton_raphson (float (*f)(float), float (*df)(float), float x0, float epsilon, int max_iter) {
    float x = x0;
    for (int i = 0; i < max_iter; i++) {
        float fx = f (x);
        float dfx = df (x);
        if (fabsf (dfx) < epsilon) {return 0.0 / 0.0;}
        float x_next = x - fx / dfx;
        if (fabsf (x_next - x) < epsilon) {return x_next;}
        x = x_next;
    } return x;
} inline float bisection (float (*f)(float), float a, float b, float epsilon, int max_iter);
inline float bisection (float (*f)(float), float a, float b, float epsilon, int max_iter) {
    float fa = f (a);
    float fb = f (b);
    if (fa * fb > 0) {return 0.0 / 0.0;}
    float c = a;
    for (int i = 0; i < max_iter; i++) {
        c = (a + b) / 2;
        float fc = f (c);
        if (fabsf (fc) < epsilon || fabsf (b - a) < epsilon) {return c;}
        if (fa * fc <= 0) {b = c; fb = fc;}
        else {a = c; fa = fc;}
    } return c;
} inline float derivative (float (*f)(float), float x, float h);
inline float derivative (float (*f)(float), float x, float h) {
    return (f (x + h) - f (x - h)) / (2 * h);
} inline float second_derivative (float (*f)(float), float x, float h);
inline float second_derivative (float (*f)(float), float x, float h) {
    return (f (x + h) - 2 * f (x) + f (x - h)) / (h * h);
} inline float integral_simpson (float (*f)(float), float a, float b, int n);
inline float integral_simpson (float (*f)(float), float a, float b, int n) {
    if (n % 2 == 1) {n++;}
    float h = (b - a) / n;
    float result = f (a) + f (b);
    for (int i = 1; i < n; i++) {
        float x = a + i * h;
        if (i % 2 == 0) {result += 2 * f (x);}
        else {result += 4 * f (x);}
    } return result * h / 3;
}

#endif
