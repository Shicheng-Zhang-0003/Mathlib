#ifndef LIBMATHC_OPTIMIZATION_H
#define LIBMATHC_OPTIMIZATION_H

#include <math.h>
#include "numerical.h"

typedef double (*opt_func)(double x);

inline double optimize_golden(opt_func f, double a, double b, double tol, int max_iter) {
    double phi = (1.0 + sqrt(5.0)) / 2.0;
    double resphi = 2.0 - phi;

    double x1 = a + resphi * (b - a);
    double x2 = b - resphi * (b - a);
    double f1 = f(x1);
    double f2 = f(x2);

    for (int i = 0; i < max_iter; i++) {
        if (b - a < tol) break;
        if (f1 < f2) {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = a + resphi * (b - a);
            f1 = f(x1);
        } else {
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = b - resphi * (b - a);
            f2 = f(x2);
        }
    }
    return (a + b) / 2.0;
}

inline double optimize_gradient_descent(opt_func f, double start, double lr, double tol, int max_iter) {
    double x = start;
    for (int i = 0; i < max_iter; i++) {
        double grad = derivative(f, x, 1e-5);
        double x_new = x - lr * grad;
        if (fabs(x_new - x) < tol) break;
        x = x_new;
    }
    return x;
}
#endif
