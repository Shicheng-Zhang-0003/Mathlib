#ifndef MATHLIB_CPU_DISPATCH_H
#define MATHLIB_CPU_DISPATCH_H

#include "ml_compiler.h"

/* ============================================================================
 * v11S CPU DISPATCH & CAPABILITIES
 *
 * ARCHITECTURAL SHIFT: Runtime dispatch and global state have been eradicated
 * to guarantee thread-safety and strict adherence to the "No Global State"
 * contract. Dispatch is now handled at compile-time via ML_COMPILE_TIME_AVX2.
 * ========================================================================== */

/* Core Matrix Multiplication API (Exported) */
ML_API void ml_matmul(const double* ML_RESTRICT A, const double* ML_RESTRICT B, double* ML_RESTRICT C, int N);

/* Compile-Time Capability Queries (Exported) */
ML_API int ml_cpu_has_avx2(void);
ML_API int ml_cpu_has_fma(void);
ML_API int ml_cpu_has_sse41(void);
ML_API int ml_cpu_has_neon(void);

#endif /* MATHLIB_CPU_DISPATCH_H */
