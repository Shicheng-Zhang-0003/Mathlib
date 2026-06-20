#ifndef LIBMATHC_QUADRATICS_H
#define LIBMATHC_QUADRATICS_H

//Library Header Files for Quadratics
#include <math.h>
inline float equation (float a, float b, float c, float x);
inline float equation (float a, float b, float c, float x) {return (a * x * x) + b * x + c;}
inline float formula_pos (float a, float b, float c);
inline float formula_pos (float a, float b, float c) {return (-b + sqrt (b * b - 4 * a * c)) / (2 * a);}
inline float formula_neg (float a, float b, float c);
inline float formula_neg (float a, float b, float c) {return (-b - sqrt (b * b - 4 * a * c)) / (2 * a);}


#endif
