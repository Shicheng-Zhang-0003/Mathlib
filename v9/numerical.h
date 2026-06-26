#ifndef LIBMATHC_NUMERICAL_H
#define LIBMATHC_NUMERICAL_H

//Library header file for numerical methods
#include <math.h>
inline double newton_raphson (double (*f)(double), double (*df)(double), double x0, double epsilon, int max_iter);
inline double newton_raphson (double (*f)(double), double (*df)(double), double x0, double epsilon, int max_iter) {
    double x = x0;
    for (int i = 0; i < max_iter; i++) {
        double fx = f (x);
        double dfx = df (x);
        if (fabs (dfx) < epsilon) {return 0.0 / 0.0;}
        double x_next = x - fx / dfx;
        if (fabs (x_next - x) < epsilon) {return x_next;}
        x = x_next;
    } return x;
} inline double bisection (double (*f)(double), double a, double b, double epsilon, int max_iter);
inline double bisection (double (*f)(double), double a, double b, double epsilon, int max_iter) {
    double fa = f (a);
    double fb = f (b);
    if (fa * fb > 0) {return 0.0 / 0.0;}
    double c = a;
    for (int i = 0; i < max_iter; i++) {
        c = (a + b) / 2;
        double fc = f (c);
        if (fabs (fc) < epsilon || fabs (b - a) < epsilon) {return c;}
        if (fa * fc <= 0) {b = c; fb = fc;}
        else {a = c; fa = fc;}
    } return c;
} inline double derivative (double (*f)(double), double x, double h);
inline double derivative (double (*f)(double), double x, double h) {
    return (f (x + h) - f (x - h)) / (2 * h);
} inline double second_derivative (double (*f)(double), double x, double h);
inline double second_derivative (double (*f)(double), double x, double h) {
    return (f (x + h) - 2 * f (x) + f (x - h)) / (h * h);
} inline double integral_simpson (double (*f)(double), double a, double b, int n);
inline double integral_simpson (double (*f)(double), double a, double b, int n) {
    if (n % 2 == 1) {n++;}
    double h = (b - a) / n;
    double result = f (a) + f (b);
    for (int i = 1; i < n; i++) {
        double x = a + i * h;
        if (i % 2 == 0) {result += 2 * f (x);}
        else {result += 4 * f (x);}
    } return result * h / 3;
}

#endif
