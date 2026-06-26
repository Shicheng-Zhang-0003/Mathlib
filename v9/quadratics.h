#ifndef LIBMATHC_QUADRATICS_H
#define LIBMATHC_QUADRATICS_H

//Library Header Files for Quadratics
#include <math.h>
inline double equation (double a, double b, double c, double x);
inline double equation (double a, double b, double c, double x) {return (a * x * x) + b * x + c;}
inline double formula_pos (double a, double b, double c);
inline double formula_pos (double a, double b, double c) {return (-b + sqrt (b * b - 4 * a * c)) / (2 * a);}
inline double formula_neg (double a, double b, double c);
inline double formula_neg (double a, double b, double c) {return (-b - sqrt (b * b - 4 * a * c)) / (2 * a);}


#endif
