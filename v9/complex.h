#ifndef LIBMATHC_COMPLEX_H
#define LIBMATHC_COMPLEX_H

//Library header file for complex numbers
#include <math.h>
#include "exponential.h"
#include "trigonometry.h"
typedef struct { double real; double imag; } cplx;
inline cplx cplx_add (cplx a, cplx b);
inline cplx cplx_add (cplx a, cplx b) {return (cplx) {a.real + b.real, a.imag + b.imag};}
inline cplx cplx_sub (cplx a, cplx b);
inline cplx cplx_sub (cplx a, cplx b) {return (cplx) {a.real - b.real, a.imag - b.imag};}
inline cplx cplx_mul (cplx a, cplx b);
inline cplx cplx_mul (cplx a, cplx b) {
    return (cplx) {a.real * b.real - a.imag * b.imag, a.real * b.imag + a.imag * b.real};
} inline cplx cplx_div (cplx a, cplx b);
inline cplx cplx_div (cplx a, cplx b) {
    double denom = b.real * b.real + b.imag * b.imag;
    if (denom == 0.0) return (cplx){0.0/0.0, 0.0/0.0};
    return (cplx) {(a.real * b.real + a.imag * b.imag) / denom, (a.imag * b.real - a.real * b.imag) / denom};
} inline double cplx_abs (cplx a);
inline double cplx_abs (cplx a) {return sqrt (a.real * a.real + a.imag * a.imag);}
inline double cplx_arg (cplx a);
inline double cplx_arg (cplx a) {
    if (a.real > 0) {return arctangent (a.imag / a.real);}
    if (a.real < 0 && a.imag >= 0) {return arctangent (a.imag / a.real) + math_pi;}
    if (a.real < 0 && a.imag < 0) {return arctangent (a.imag / a.real) - math_pi;}
    if (a.real == 0 && a.imag > 0) {return math_pi / 2;}
    if (a.real == 0 && a.imag < 0) {return -math_pi / 2;}
    return 0.0 / 0.0;
} inline cplx cplx_conjugate (cplx a);
inline cplx cplx_conjugate (cplx a) {return (cplx) {a.real, -a.imag};}
inline cplx cplx_exponential (cplx a);
inline cplx cplx_exponential (cplx a) {
    double mag = exponential (a.real);
    return (cplx) {mag * cosine (a.imag), mag * sine (a.imag)};
} inline cplx cplx_logarithm (cplx a);
inline cplx cplx_logarithm (cplx a) {
    return (cplx) {logarithm (cplx_abs (a)), cplx_arg (a)};
} inline cplx cplx_power (cplx a, cplx b);
inline cplx cplx_power (cplx a, cplx b) {
    cplx log_a = cplx_logarithm (a);
    return cplx_exponential ((cplx) {b.real * log_a.real - b.imag * log_a.imag, b.real * log_a.imag + b.imag * log_a.real});
}



#endif
