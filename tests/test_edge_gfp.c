#include "test_harness.h"
#include "ml_gfp.h"

#include <limits.h>
#include <stdint.h>
int main (void) {
    ml_test_ctx_t ctx;
    ml_gfp_ctx_t f7, f97, f2, large;
    uint32_t r = 123u, check;
    int symbol = 7;
    ml_test_init (&ctx, "Prime finite-field foundation");

    ml_gfp_ctx_init (&f7);
    ml_gfp_ctx_init (&f97);
    ml_gfp_ctx_init (&f2);
    ml_gfp_ctx_init (&large);
    ASSERT_TRUE (&ctx, ml_gfp_ctx_set_prime (&f7, 7u) == ML_SUCCESS && ml_gfp_ctx_is_valid (&f7), "construct prime field F7");
    ASSERT_TRUE (&ctx, ml_gfp_ctx_set_prime (&f7, 91u) == ML_ERR_INVALID_ARG && f7.p == 7u && ml_gfp_ctx_is_valid (&f7),
                 "composite modulus rejected without replacing valid context");
    ASSERT_TRUE (&ctx, ml_gfp_add (&f7, 5u, 6u, &r) == ML_SUCCESS && r == 4u,
                 "field addition wraps modulo p");
    ASSERT_TRUE (&ctx, ml_gfp_sub (&f7, 1u, 6u, &r) == ML_SUCCESS && r == 2u,
                 "field subtraction wraps modulo p");
    ASSERT_TRUE (&ctx, ml_gfp_mul (&f7, 6u, 6u, &r) == ML_SUCCESS && r == 1u,
                 "field multiplication uses exact wide intermediate");
    ASSERT_TRUE (&ctx, ml_gfp_inv (&f7, 3u, &r) == ML_SUCCESS && r == 5u,
                 "multiplicative inverse by Fermat exponent");
    ASSERT_TRUE (&ctx, ml_gfp_div (&f7, 4u, 2u, &r) == ML_SUCCESS && r == 2u,
                 "field division");
    ASSERT_TRUE (&ctx, ml_gfp_pow (&f7, 0u, 0u, &r) == ML_SUCCESS && r == 1u,
                 "zero exponent convention");
    ASSERT_TRUE (&ctx, ml_gfp_from_i64 (&f7, INT64_MIN, &r) == ML_SUCCESS && r == 6u,
                 "signed normalization handles INT64_MIN");
    ASSERT_TRUE (&ctx, ml_gfp_legendre (&f7, 2u, &symbol) == ML_SUCCESS && symbol == 1,
                 "Euler criterion identifies quadratic residue");
    ASSERT_TRUE (&ctx, ml_gfp_legendre (&f7, 3u, &symbol) == ML_SUCCESS && symbol == -1,
                 "Euler criterion identifies nonsquare");
    ASSERT_TRUE (&ctx, ml_gfp_legendre (&f7, 0u, &symbol) == ML_SUCCESS && symbol == 0,
                 "Legendre symbol of zero");
    ASSERT_TRUE (&ctx, ml_gfp_sqrt (&f7, 2u, &r) == ML_SUCCESS && r == 3u,
                 "square root in the p equals 3 mod 4 branch");
    r = 111u;
    ASSERT_TRUE (&ctx, ml_gfp_sqrt (&f7, 3u, &r) == ML_ERR_SINGULAR && r == 111u,
                 "nonsquare root fails without modifying output");
    r = 111u;
    ASSERT_TRUE (&ctx, ml_gfp_div (&f7, 1u, 0u, &r) == ML_ERR_SINGULAR && r == 111u,
                 "division by zero fails without modifying output");
    ASSERT_TRUE (&ctx, ml_gfp_add (&f7, 7u, 0u, &r) == ML_ERR_INVALID_ARG && r == 111u,
                 "noncanonical residue rejected atomically");

    ASSERT_TRUE (&ctx, ml_gfp_ctx_set_prime (&f97, 97u) == ML_SUCCESS && ml_gfp_sqrt (&f97, 12u, &r) == ML_SUCCESS && ml_gfp_mul (&f97, r, r, &check) == ML_SUCCESS && check == 12u && r == 20u,
                 "Tonelli-Shanks square root for p equals 1 mod 4");
    ASSERT_TRUE (&ctx, ml_gfp_ctx_set_prime (&f2, 2u) == ML_SUCCESS && ml_gfp_sqrt (&f2, 1u, &r) == ML_SUCCESS && r == 1u,
                 "characteristic two prime field");
    ASSERT_TRUE (&ctx, ml_gfp_ctx_set_prime (&large, UINT32_C (4294967291)) == ML_SUCCESS && ml_gfp_mul (&large, UINT32_C (4294967290), UINT32_C (4294967290), &r) == ML_SUCCESS && r == 1u,
                 "largest 32-bit prime modulus and wide product");

    return ml_test_summary (&ctx);
}
