#ifndef LIBMATHC_ERROR_FREE_H
#define LIBMATHC_ERROR_FREE_H

// Dekker's Fast Two-Sum (Assumes |a| >= |b|)
inline double ml_fast_two_sum(double a, double b, double *err) {
    double s = a + b;
    double z = s - a;
    *err = b - z;
    return s;
}

// Knuth's Two-Sum (No magnitude assumption)
inline double ml_two_sum(double a, double b, double *err) {
    double s = a + b;
    double v = s - a;
    *err = (a - (s - v)) + (b - v);
    return s;
}

// Dekker's Two-Product
inline double ml_two_product(double a, double b, double *err) {
    double p = a * b;
    // Split a and b into high/low 26-bit parts
    double ca = a * 67108865.0; // 2^26 + 1
    double a_hi = ca - (ca - a);
    double a_lo = a - a_hi;
    double cb = b * 67108865.0;
    double b_hi = cb - (cb - b);
    double b_lo = b - b_hi;
    *err = ((a_hi * b_hi - p) + a_hi * b_lo + a_lo * b_hi) + a_lo * b_lo;
    return p;
}

// Software FMA (Fused Multiply-Add) using Error-Free transformations
inline double ml_fma(double a, double b, double c) {
    double p, err;
    p = ml_two_product(a, b, &err);
    double s1, s2;
    s1 = ml_two_sum(p, c, &s2);
    return s1 + (err + s2); // Captures all intermediate rounding errors
}
#endif
