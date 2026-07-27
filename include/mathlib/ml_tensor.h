#ifndef MATHLIB_ML_TENSOR_H
#define MATHLIB_ML_TENSOR_H

#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include "ml_compiler.h"
#include "ml_types.h"

/* ============================================================================
 * v11S TENSOR & WORKSPACE SYSTEM
 *
 * SAFETY BY DEFAULT: All bounds checks, NULL checks, and dimension validation
 * are UNCONDITIONAL. The MATHLIB_PROFILE_HARDENED flag now ONLY controls
 * memory canaries and NaN poisoning, NOT basic safety.
 *
 * ABI STABILITY: The ml_workspace_t struct has a FIXED size regardless of
 * compile flags. The magic_canary is ALWAYS present to guarantee binary
 * compatibility across translation units.
 * ========================================================================== */

typedef struct {
    uint8_t *base;
    size_t capacity;
    size_t offset;
    uint64_t magic_canary;
} ml_workspace_t;

#define ML_WORKSPACE_CANARY 0xDEADBEEFCAFEBABEULL

ML_INLINE void ml_workspace_init(ml_workspace_t *ws, void *buffer, size_t size) {
    if (ML_UNLIKELY(!ws)) return;

    if (ML_UNLIKELY(!buffer || size == 0)) {
        ws->base = NULL;
        ws->capacity = 0;
        ws->offset = 0;
        ws->magic_canary = 0;
        return;
    }

    uintptr_t raw = (uintptr_t)buffer;
    uintptr_t aligned = (raw + (uintptr_t)31) & ~(uintptr_t)31;
    size_t pad = (size_t)(aligned - raw);

    if (ML_UNLIKELY(pad >= size)) {
        ws->base = NULL;
        ws->capacity = 0;
        ws->offset = 0;
        ws->magic_canary = 0;
        return;
    }

    ws->base = (uint8_t *)aligned;
    ws->capacity = size - pad;
    ws->offset = 0;
    ws->magic_canary = ML_WORKSPACE_CANARY;
}

ML_INLINE void *ml_workspace_alloc(ml_workspace_t *ws, size_t bytes) {
    if (ML_UNLIKELY(!ws)) return NULL;
    if (ML_UNLIKELY(ws->magic_canary != ML_WORKSPACE_CANARY)) return NULL;

    size_t aligned = (bytes + 31) & ~(size_t)31;

    if (ML_UNLIKELY(aligned < bytes)) return NULL;
    if (ML_UNLIKELY(ws->offset > ws->capacity)) return NULL;
    if (ML_UNLIKELY(aligned > ws->capacity - ws->offset)) return NULL;

    void *ptr = ws->base + ws->offset;
    ws->offset += aligned;
    return ptr;
}

ML_INLINE void ml_workspace_reset(ml_workspace_t *ws) {
    if (ML_UNLIKELY(!ws)) return;
    if (ML_UNLIKELY(ws->magic_canary != ML_WORKSPACE_CANARY)) return;
    ws->offset = 0;
}

typedef struct {
    double *data;
    int rows;
    int cols;
} ml_tensor_view_t;

ML_INLINE ml_tensor_view_t ml_tensor_view(double *data, int rows, int cols) {
    return (ml_tensor_view_t){data, rows, cols};
}

/* MATHLIB_CLOSURE_P1_SIZE_T_INDEXING */
ML_INLINE double *ml_tensor_at(ml_tensor_view_t t, int r, int c) {
    if (ML_UNLIKELY(!t.data)) return NULL;
    if (ML_UNLIKELY(r < 0 || r >= t.rows || c < 0 || c >= t.cols)) return NULL;

    return &t.data[(size_t)r * (size_t)t.cols + (size_t)c];
}

#define ML_TENSOR_AT(t, r, c) ((t).data[(size_t)(r) * (size_t)(t).cols + (size_t)(c)])

#endif /* MATHLIB_ML_TENSOR_H */
