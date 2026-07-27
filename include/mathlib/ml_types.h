#ifndef MATHLIB_TYPES_H
#define MATHLIB_TYPES_H
#include "ml_core.h"

#include <stdint.h>
#include <stddef.h>

/* ============================================================================
 * v11S STRICT TYPE SYSTEM
 *
 * Public API continues to use standard C types (double, float, int) for
 * maximum interoperability. Internal implementations and hardware-specific
 * kernels should use these strict typedefs to guarantee exact bit-widths
 * across all platforms (especially for the Yellow Team's MSVC targets).
 * ========================================================================== */

/* Floating Point Types */
typedef float  ml_f32;
typedef double ml_f64;

/* Signed Integer Types */
typedef int8_t   ml_i8;
typedef int16_t  ml_i16;
typedef int32_t  ml_i32;
typedef int64_t  ml_i64;

/* Unsigned Integer Types */
typedef uint8_t  ml_u8;
typedef uint16_t ml_u16;
typedef uint32_t ml_u32;
typedef uint64_t ml_u64;

/* Pointer & Size Types */
typedef size_t    ml_size;
typedef ptrdiff_t ml_ptrdiff;

/* ============================================================================
 * v11S UNIFIED STATUS CODES
 * ========================================================================== */
typedef enum {
    ML_SUCCESS          =  0,
    ML_ERR_SINGULAR     = -1,
    ML_ERR_WORKSPACE    = -2,
    ML_ERR_INVALID_ARG  = -3,
    ML_ERR_NAN_INPUT    = -4,
    ML_ERR_INTERNAL     = -99
} ml_status_t;

#endif /* MATHLIB_TYPES_H */
