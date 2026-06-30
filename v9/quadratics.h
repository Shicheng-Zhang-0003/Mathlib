#ifndef LIBMATHC_QUADRATICS_H
#define LIBMATHC_QUADRATICS_H

//Library Header Files for Quadratics
#include "ml_core.h"
static inline double equation (double a, double b, double c, double x);
static inline double equation (double a, double b, double c, double x) {return (a * x * x) + b * x + c;}
static inline double formula_pos (double a, double b, double c);
static inline double formula_pos (double a, double b, double c) {return (-b + ml_sqrt(b * b - 4 * a * c)) / (2 * a);}
static inline double formula_neg (double a, double b, double c);
static inline double formula_neg (double a, double b, double c) {return (-b - ml_sqrt(b * b - 4 * a * c)) / (2 * a);}


#endif
