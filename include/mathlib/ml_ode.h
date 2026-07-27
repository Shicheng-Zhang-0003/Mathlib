#ifndef MATHLIB_ML_ODE_H
#define MATHLIB_ML_ODE_H

#include "ml_compiler.h"
#include "ml_core.h"

typedef double (*ml_ode_func_t)(double t, double y);

ML_API double ml_ode_euler(ml_ode_func_t f, double t0, double y0, double dt, int steps);
ML_API double ml_ode_rk4(ml_ode_func_t f, double t0, double y0, double dt, int steps);

#endif /* MATHLIB_ML_ODE_H */
