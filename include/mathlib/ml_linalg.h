#ifndef MATHLIB_V10_LINALG_H
#define MATHLIB_V10_LINALG_H

#include "ml_compiler.h"
#include "ml_tensor.h"
#include "ml_core.h"
#include "ml_types.h"

ML_API ml_status_t ml_lu_decomp(ml_tensor_view_t A, ml_tensor_view_t LU, int* P, ml_workspace_t* ws);
ML_API ml_status_t ml_solve(ml_tensor_view_t A, double* b, double* x, ml_workspace_t* ws);


/* Matrix-Vector Multiplication (y = Ax) */
ML_API void ml_matvec(ml_tensor_view_t A, const double* x, double* out);

#endif /* MATHLIB_V10_LINALG_H */
