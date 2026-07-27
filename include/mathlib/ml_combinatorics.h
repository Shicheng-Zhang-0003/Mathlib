#ifndef MATHLIB_ML_COMBINATORICS_H
#define MATHLIB_ML_COMBINATORICS_H

#include <stdint.h>
#include "ml_compiler.h"

ML_API uint64_t ml_factorial(int x);
ML_API uint64_t ml_npr(int n, int r);
ML_API uint64_t ml_ncr(int n, int r);

#endif /* MATHLIB_ML_COMBINATORICS_H */
