#ifndef MATHLIB_TEST_HARNESS_H
#define MATHLIB_TEST_HARNESS_H

#include <stdio.h>
#include <stdlib.h>
#include "ml_core.h"

/* v11S Test Context: Eradicates static globals to prevent linker bombs
 * and allow multi-suite aggregation if needed. */
typedef struct {
    int passed;
    int failed;
    const char* suite_name;
} ml_test_ctx_t;

static inline void ml_test_init(ml_test_ctx_t* ctx, const char* name) {
    ctx->passed = 0;
    ctx->failed = 0;
    ctx->suite_name = name;
}

#define ASSERT_TRUE(ctx, cond, msg) do { \
    if (cond) { (ctx)->passed++; } \
    else { (ctx)->failed++; printf("  [FAIL] %s (Line %d)\n", msg, __LINE__); } \
} while(0)

#define ASSERT_NEAR(ctx, a, b, eps, msg) do { \
    double _a = (double)(a); double _b = (double)(b); \
    if (ml_fabs(_a - _b) <= (eps)) { (ctx)->passed++; } \
    else { (ctx)->failed++; printf("  [FAIL] %s (Line %d) got %.17g expected %.17g diff %.17g\n", msg, __LINE__, _a, _b, ml_fabs(_a - _b)); } \
} while(0)

static inline int ml_test_summary(ml_test_ctx_t* ctx) {
    printf("[%s] Passed: %d, Failed: %d\n", ctx->suite_name, ctx->passed, ctx->failed);
    return ctx->failed > 0 ? 1 : 0;
}

#endif /* MATHLIB_TEST_HARNESS_H */
