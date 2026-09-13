#ifndef MATHLIB_COMPILER_H
#define MATHLIB_COMPILER_H

/* ============================================================================
 * MATHLIB v11S COMPILER ABSTRACTION LAYER
 * All compiler-specific extensions, intrinsics, and attributes must be routed
 * through these macros to ensure cross-platform compatibility (GCC/Clang/MSVC).
 * ========================================================================== */

/* 1. Compiler Identity */
#if defined(_MSC_VER)
#  define ML_COMPILER_MSVC 1
#elif defined(__clang__)
#  define ML_COMPILER_CLANG 1
#elif defined(__GNUC__)
#  define ML_COMPILER_GNUC 1
#else
#  define ML_COMPILER_UNKNOWN 1
#endif

/* 2. Inline & Alignment Abstraction */
#if defined(ML_COMPILER_MSVC)
#  define ML_INLINE static __forceinline
#  define ML_NOINLINE __declspec(noinline)
#  define ML_ALIGN(n) __declspec(align(n))
#else
#  define ML_INLINE static inline __attribute__((always_inline))
#  define ML_NOINLINE __attribute__((noinline))
#  define ML_ALIGN(n) __attribute__((aligned(n)))
#endif

/* 3. API Export (Harmless for static libs, required for future DLLs) */
#if defined(ML_COMPILER_MSVC)
#  if defined(MATHLIB_BUILD_SHARED)
#    define ML_API __declspec(dllexport)
#  elif defined(MATHLIB_USE_SHARED)
#    define ML_API __declspec(dllimport)
#  else
#    define ML_API
#  endif
#elif defined(ML_COMPILER_GNUC) || defined(ML_COMPILER_CLANG)
#  define ML_API __attribute__((visibility("default")))
#else
#  define ML_API
#endif

/* 4. Branch Prediction Hints */
#if defined(ML_COMPILER_GNUC) || defined(ML_COMPILER_CLANG)
#  define ML_LIKELY(x) __builtin_expect(!!(x), 1)
#  define ML_UNLIKELY(x) __builtin_expect(!!(x), 0)
#else
#  define ML_LIKELY(x) (x)
#  define ML_UNLIKELY(x) (x)
#endif

/* 5. Restrict Pointer (C99 / MSVC) */
#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)
#  define ML_RESTRICT restrict
#elif defined(ML_COMPILER_MSVC)
#  define ML_RESTRICT __restrict
#else
#  define ML_RESTRICT
#endif

/* 6. Target Attributes (for AVX2/FMA routing) */
#if defined(ML_COMPILER_GNUC) || defined(ML_COMPILER_CLANG)
#  define ML_TARGET_AVX2 __attribute__((target("avx2,fma")))
#else
#  define ML_TARGET_AVX2 /* MSVC uses /arch:AVX2 globally or pragmas */
#endif

/* 7. Constructor (Library auto-init) */
#if defined(ML_COMPILER_GNUC) || defined(ML_COMPILER_CLANG)
#  define ML_CTOR __attribute__((constructor))
#else
#  define ML_CTOR /* MSVC requires #pragma init_seg or DllMain */
#endif

/* 8. Deprecation Warnings */
#if defined(ML_COMPILER_GNUC) || defined(ML_COMPILER_CLANG)
#  define ML_DEPRECATED(msg) __attribute__((deprecated(msg)))
#elif defined(ML_COMPILER_MSVC)
#  define ML_DEPRECATED(msg) __declspec(deprecated(msg))
#else
#  define ML_DEPRECATED(msg)
#endif

/* 9. Compile-time Feature Detection */
#if defined(__AVX2__) && (defined(__x86_64__) || defined(__i386__))
#  define ML_COMPILE_TIME_AVX2 1
#else
#  define ML_COMPILE_TIME_AVX2 0
#endif


/* 10. Hardware FMA (Fused Multiply-Add) Abstraction */
#if defined(__GNUC__) || defined(__clang__)
#  define ML_FMA(a, b, c) __builtin_fma((double)(a), (double)(b), (double)(c))
#elif defined(_MSC_VER)
#  include <math.h>
#  define ML_FMA(a, b, c) fma((double)(a), (double)(b), (double)(c))
#else
/* Software FMA emulation using Dekker's Two-Product + Two-Sum.
 * This provides single-rounding semantics (approximately) on platforms without FMA.
 * For production use, hardware FMA is strongly recommended. */
static inline double ml_fma_soft_impl(double a, double b, double c) {
    double p = a * b;
    double c2 = c;
    /* Two-Product: p = fl(a*b), err = a*b - p */
    double ca = a * 134217729.0;
    double a_hi = ca - (ca - a);
    double a_lo = a - a_hi;
    double cb = b * 134217729.0;
    double b_hi = cb - (cb - b);
    double b_lo = b - b_hi;
    double err = ((a_hi * b_hi - p) + a_hi * b_lo + a_lo * b_hi) + a_lo * b_lo;
    /* Two-Sum: s = fl(p + c), err2 = p + c - s */
    double s = p + c2;
    double v = s - p;
    double err2 = (p - (s - v)) + (c2 - v);
    return s + (err + err2);
}
#  define ML_FMA(a, b, c) ml_fma_soft_impl((double)(a), (double)(b), (double)(c))
#endif

#endif /* MATHLIB_COMPILER_H */
