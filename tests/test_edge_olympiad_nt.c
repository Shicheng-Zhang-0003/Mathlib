/* Olympiad-NT + recreational pins: LTE, floor-sum, Legendre, CF, Pell,
 * digits, Collatz, magic squares. Ground truth: hand/Python exact. */
#include "test_harness.h"
#include "ml_numbertheory.h"
#include "ml_contest.h"
#include "ml_core.h"
#include <stdint.h>
int main (void) {
    ml_test_ctx_t ctx;
    ml_test_init (&ctx, "Edge Olympiad NT");

    /* Legendre n!. */
    ASSERT_TRUE (&ctx, ml_vp_factorial (100, 5) == 24, "v5(100!)");
    ASSERT_TRUE (&ctx, ml_vp_factorial (10, 2) == 8, "v2(10!)");
    ASSERT_TRUE (&ctx, ml_vp_factorial (0, 5) == 0, "v5(0!)");
    ASSERT_TRUE (&ctx, ml_vp_factorial (100, 6) == UINT64_MAX, "non-prime sentinel");
    ASSERT_TRUE (&ctx, ml_vp_factorial (-3, 5) == UINT64_MAX, "negative sentinel");

    /* Floor sum. */
    ASSERT_TRUE (&ctx, ml_floor_sum_u64 (10, 7, 3, 2) == 18, "floor_sum basic");
    ASSERT_TRUE (&ctx, ml_floor_sum_u64 (0, 7, 3, 2) == 0, "floor_sum n=0");
    ASSERT_TRUE (&ctx, ml_floor_sum_u64 (10, 0, 3, 2) == UINT64_MAX, "floor_sum m=0");

    /* LTE. */
    ASSERT_TRUE (&ctx, ml_lte_vp (10, 4, 6, 3) == 2, "LTE 3-adic");
    ASSERT_TRUE (&ctx, ml_lte_vp (15, 8, 3, 7) == 1, "LTE 7-adic");
    ASSERT_TRUE (&ctx, ml_lte_vp (5, 3, 4, 2) == 5, "LTE 2-adic even");
    ASSERT_TRUE (&ctx, ml_lte_vp (5, 3, 3, 2) == 1, "LTE 2-adic odd (v2(125-27)=1)");
    ASSERT_TRUE (&ctx, ml_lte_vp (10, 4, 6, 5) == UINT64_MAX, "LTE p not dividing");
    ASSERT_TRUE (&ctx, ml_lte_vp (5, 5, 3, 7) == UINT64_MAX, "LTE a==b");

    /* Continued fractions + Pell. */
    ASSERT_TRUE (&ctx, ml_sqrt_cf_period (2) == 1, "cf sqrt2");
    ASSERT_TRUE (&ctx, ml_sqrt_cf_period (23) == 4, "cf sqrt23");
    ASSERT_TRUE (&ctx, ml_sqrt_cf_period (61) == 11, "cf sqrt61");
    ASSERT_TRUE (&ctx, ml_sqrt_cf_period (4) == 0, "cf square");
    ASSERT_TRUE (&ctx, ml_sqrt_cf_period (1) == -1, "cf invalid");
    {
        uint64_t x = 0, y = 0, n = 0, d = 0;
        ASSERT_TRUE (&ctx, ml_pell_fund (2, &x, &y) == ML_SUCCESS && x == 3 && y == 2, "Pell D=2");
        ASSERT_TRUE (&ctx, ml_pell_fund (5, &x, &y) == ML_SUCCESS && x == 9 && y == 4, "Pell D=5");
        ASSERT_TRUE (&ctx, ml_pell_fund (61, &x, &y) == ML_SUCCESS, "Pell D=61 status");
        ASSERT_TRUE (&ctx, x == 1766319049ULL && y == 226153980ULL, "Pell D=61 famous");
        ASSERT_TRUE (&ctx, ml_pell_fund (4, &x, &y) == ML_ERR_INVALID_ARG, "Pell square");
        ASSERT_TRUE (&ctx, ml_cf_convergent_sqrt (2, 0, &n, &d) == ML_SUCCESS && n == 1 && d == 1, "conv0");
        ASSERT_TRUE (&ctx, ml_cf_convergent_sqrt (2, 3, &n, &d) == ML_SUCCESS && n == 17 && d == 12, "conv3");
    }

    /* Recreational. */
    ASSERT_TRUE (&ctx, ml_digit_sum_u64 (12345) == 15, "digit sum");
    ASSERT_TRUE (&ctx, ml_digit_sum_u64 (0) == 0, "digit sum 0");
    ASSERT_TRUE (&ctx, ml_digital_root (12345) == 6, "digital root");
    ASSERT_TRUE (&ctx, ml_digital_root (9) == 9, "digital root 9");
    ASSERT_TRUE (&ctx, ml_is_palindrome_u64 (12321) == 1, "palindrome");
    ASSERT_TRUE (&ctx, ml_is_palindrome_u64 (12345) == 0, "non-palindrome");
    ASSERT_TRUE (&ctx, ml_collatz_steps (6) == 8, "collatz 6");
    ASSERT_TRUE (&ctx, ml_collatz_steps (1) == 0, "collatz 1");
    ASSERT_TRUE (&ctx, ml_collatz_steps (27) == 111, "collatz 27");
    ASSERT_TRUE (&ctx, ml_collatz_steps (0) == UINT64_MAX, "collatz 0");
    {
        long long lo [] = {8, 1, 6, 3, 5, 7, 4, 9, 2};
        long long bad [] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
        ASSERT_TRUE (&ctx, ml_magic_square_check (lo, 3) == 1, "Lo Shu");
        ASSERT_TRUE (&ctx, ml_magic_square_check (bad, 3) == 0, "non-magic");
        ASSERT_TRUE (&ctx, ml_magic_square_check (lo, 2) == 0, "magic too small");
    }

    return ml_test_summary (&ctx);
}
