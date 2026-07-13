#ifndef MATHLIB_V10_TENSOR_H
#define MATHLIB_V10_TENSOR_H

#include <stdint.h>
#include <stddef.h>

// A lightweight view into externally managed memory.
// No ownership, no allocations. Allows strided slicing.
typedef struct {
    double* data;
    int rows;
    int cols;
    int row_stride; // Distance in elements between rows (allows sub-matrix views)
} ml_tensor_view_t;

// A pre-allocated scratchpad for solvers.
// Replaces internal malloc/free calls in hot loops.
typedef struct {
    void* storage;
    size_t size_bytes;
    size_t used_bytes;
#if defined(MATHLIB_PROFILE_HARDENED)
    uint64_t magic_canary;
#endif
} ml_workspace_t;

// Bump allocator for the workspace (32-byte aligned for AVX)
static inline void ml_workspace_init(ml_workspace_t* ws) {
    ws->used_bytes = 0;
#if defined(MATHLIB_PROFILE_HARDENED)
    ws->magic_canary = 0xDEADBEEFCAFEBABEULL;
#endif
}

static inline void* ml_workspace_alloc(ml_workspace_t* ws, size_t bytes) {
#if defined(MATHLIB_PROFILE_HARDENED)
    if (ws->magic_canary != 0xDEADBEEFCAFEBABEULL) return NULL; // Structural corruption detected
#endif
    size_t aligned = (bytes + 31) & ~(size_t)31;
    if (ws->used_bytes + aligned > ws->size_bytes) return NULL; // Strict bounds assert
    void* ptr = (char*)ws->storage + ws->used_bytes;
    ws->used_bytes += aligned;
    return ptr;
}

static inline void ml_workspace_reset(ml_workspace_t* ws) {
    ws->used_bytes = 0;
}

// Helper to create a view from a raw contiguous buffer
static inline ml_tensor_view_t ml_tensor_view(double* data, int rows, int cols) {
    return (ml_tensor_view_t){data, rows, cols, cols};
}

// Helper to access elements respecting stride
#define ML_TENSOR_AT(t, r, c) ((t).data[(r) * (t).row_stride + (c)])

#endif
