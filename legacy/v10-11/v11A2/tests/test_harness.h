#ifndef MATHLIB_TEST_HARNESS_H
#define MATHLIB_TEST_HARNESS_H
#include <stdio.h>
#include "ml_core.h"
#include <stdlib.h>

static int g_passed = 0;
static int g_failed = 0;

#define ASSERT_TRUE(cond, msg) do { \
    if (cond) { g_passed++; } \
    else { g_failed++; printf("  [FAIL] %s (Line %d)\n", msg, __LINE__); } \
} while(0)

#define ASSERT_NEAR(a, b, eps, msg) ASSERT_TRUE(ml_fabs((double)(a) - (double)(b)) < (eps), msg)

static inline int test_harness_summary(const char* suite_name) {
    printf("[%s] Passed: %d, Failed: %d\n", suite_name, g_passed, g_failed);
    return g_failed > 0 ? 1 : 0;
}
#endif