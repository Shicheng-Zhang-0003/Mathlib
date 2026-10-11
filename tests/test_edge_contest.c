/* Contest extravaganza pins: figurate -> olympiad -> Putnam.
 * Ground truth: Python big-int exact (u64), mpmath 80-dps (doubles).
 * u64 paths are exact-or-sentinel; doubles are exactness-guarded or
 * LD regression-guard tier (see ml_contest.h + PRECISION_CONTRACT rule 3). */
#include "test_harness.h"
#include "ml_contest.h"
#include "ml_core.h"
#include <stdint.h>

int main(void) {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Edge Contest");

    /* Tier 1: figurate + parity + recurrences + games. */
    ASSERT_TRUE(&ctx, ml_triangular(0) == 0, "T0");
    ASSERT_TRUE(&ctx, ml_triangular(10) == 55, "T10");
    ASSERT_TRUE(&ctx, ml_triangular(100) == 5050, "T100");
    ASSERT_TRUE(&ctx, ml_pentagonal(5) == 35, "P5");
    ASSERT_TRUE(&ctx, ml_pentagonal(10) == 145, "P10");
    ASSERT_TRUE(&ctx, ml_hexagonal(5) == 45, "H5");
    ASSERT_TRUE(&ctx, ml_hexagonal(10) == 190, "H10");
    ASSERT_TRUE(&ctx, ml_is_triangular(55) == 1, "55 triangular");
    ASSERT_TRUE(&ctx, ml_is_triangular(56) == 0, "56 not triangular");
    ASSERT_TRUE(&ctx, ml_is_triangular(0) == 1, "0 triangular");
    ASSERT_TRUE(&ctx, ml_double_factorial(10) == 3840, "10!!");
    ASSERT_TRUE(&ctx, ml_double_factorial(11) == 10395, "11!!");
    ASSERT_TRUE(&ctx, ml_double_factorial(0) == 1, "0!!");
    ASSERT_TRUE(&ctx, ml_binom_odd(5, 2) == 0, "C5,2 even");
    ASSERT_TRUE(&ctx, ml_binom_odd(5, 1) == 1, "C5,1 odd");
    ASSERT_TRUE(&ctx, ml_binom_odd(7, 3) == 1, "C7,3=35 odd");
    ASSERT_TRUE(&ctx, ml_stars_bars(3, 2) == 10, "stars-bars 3,2");
    ASSERT_TRUE(&ctx, ml_hanoi_moves(0) == 0, "hanoi 0");
    ASSERT_TRUE(&ctx, ml_hanoi_moves(10) == 1023, "hanoi 10");
    ASSERT_TRUE(&ctx, ml_hanoi_moves(64) == UINT64_MAX, "hanoi 64 sentinel");
    ASSERT_TRUE(&ctx, ml_josephus(7, 3) == 3, "josephus 7,3 (0-idx)");
    ASSERT_TRUE(&ctx, ml_josephus(1, 1) == 0, "josephus 1,1");
    ASSERT_TRUE(&ctx, ml_lucas_number(0) == 2, "L0");
    ASSERT_TRUE(&ctx, ml_lucas_number(10) == 123, "L10");
    ASSERT_TRUE(&ctx, ml_lucas_number(11) == 199, "L11");
    ASSERT_TRUE(&ctx, ml_pell_number(7) == 169, "Pell7");
    ASSERT_TRUE(&ctx, ml_pell_number(10) == 2378, "Pell10");
    ASSERT_TRUE(&ctx, ml_tribonacci(10) == 81, "Trib10");
    ASSERT_TRUE(&ctx, ml_tribonacci(11) == 149, "Trib11");

    /* Tier 2: Bell/Stirling/Eulerian/Lah/Motzkin/Schroder/Narayana. */
    ASSERT_TRUE(&ctx, ml_bell_number(0) == 1, "B0");
    ASSERT_TRUE(&ctx, ml_bell_number(5) == 52, "B5");
    ASSERT_TRUE(&ctx, ml_bell_number(10) == 115975, "B10");
    ASSERT_TRUE(&ctx, ml_bell_number(12) == 4213597, "B12");
    ASSERT_TRUE(&ctx, ml_stirling2(7, 3) == 301, "S2 7,3");
    ASSERT_TRUE(&ctx, ml_stirling2(10, 5) == 42525, "S2 10,5");
    ASSERT_TRUE(&ctx, ml_stirling2(0, 0) == 1, "S2 0,0");
    ASSERT_TRUE(&ctx, ml_stirling1(7, 2) == 1764, "S1 7,2");
    ASSERT_TRUE(&ctx, ml_stirling1(7, 5) == 175, "S1 7,5");
    ASSERT_TRUE(&ctx, ml_eulerian(7, 1) == 120, "Eul 7,1");
    ASSERT_TRUE(&ctx, ml_eulerian(7, 3) == 2416, "Eul 7,3");
    ASSERT_TRUE(&ctx, ml_lah_number(6, 3) == 1200, "Lah 6,3");
    ASSERT_TRUE(&ctx, ml_motzkin(5) == 21, "Motzkin 5");
    ASSERT_TRUE(&ctx, ml_motzkin(10) == 2188, "Motzkin 10");
    ASSERT_TRUE(&ctx, ml_schroder_large(5) == 394, "Schroder 5");
    ASSERT_TRUE(&ctx, ml_schroder_large(8) == 41586, "Schroder 8");
    ASSERT_TRUE(&ctx, ml_narayana(5, 3) == 20, "Narayana 5,3");
    ASSERT_TRUE(&ctx, ml_narayana(5, 1) == 1, "Narayana 5,1");
    /* Narayana rows sum to Catalan. */
    ASSERT_TRUE(&ctx, ml_narayana(5, 1) + ml_narayana(5, 2) + ml_narayana(5, 3)
        + ml_narayana(5, 4) + ml_narayana(5, 5) == 42, "Narayana row = Catalan5");
    ASSERT_TRUE(&ctx, ml_part_distinct(5) == 3, "pdist 5");
    ASSERT_TRUE(&ctx, ml_part_distinct(10) == 10, "pdist 10");
    ASSERT_TRUE(&ctx, ml_part_distinct(20) == 64, "pdist 20");
    ASSERT_TRUE(&ctx, ml_part_distinct(50) == 3658, "pdist 50");

    /* Tier 2: prime/game/congruence staples. */
    ASSERT_TRUE(&ctx, ml_goldbach_partitions(10) == 2, "GB 10");
    ASSERT_TRUE(&ctx, ml_goldbach_partitions(28) == 2, "GB 28");
    ASSERT_TRUE(&ctx, ml_goldbach_partitions(100) == 6, "GB 100");
    ASSERT_TRUE(&ctx, ml_goldbach_partitions(9) == 0, "GB odd 0");
    ASSERT_TRUE(&ctx, ml_is_sophie_germain(11) == 1, "SG 11");
    ASSERT_TRUE(&ctx, ml_is_sophie_germain(13) == 0, "SG 13 no");
    ASSERT_TRUE(&ctx, ml_lucas_lehmer(3) == 1, "LL M3 prime");
    ASSERT_TRUE(&ctx, ml_lucas_lehmer(5) == 1, "LL M5 prime");
    ASSERT_TRUE(&ctx, ml_lucas_lehmer(11) == 0, "LL M11 composite");
    ASSERT_TRUE(&ctx, ml_lucas_lehmer(4) == -1, "LL invalid");
    ASSERT_TRUE(&ctx, ml_modsqrt_tonelli(10, 13) == 6, "sqrt10 mod13 = 6");
    ASSERT_TRUE(&ctx, ml_modsqrt_tonelli(5, 11) == 4, "sqrt5 mod11 = 4");
    ASSERT_TRUE(&ctx, ml_modsqrt_tonelli(2, 7) == 3, "sqrt2 mod7 = 3");
    ASSERT_TRUE(&ctx, ml_modsqrt_tonelli(3, 7) == UINT64_MAX, "3 nonresidue mod7");
    ASSERT_TRUE(&ctx, ml_tau_divisors(36) == 9, "tau 36");
    ASSERT_TRUE(&ctx, ml_tau_divisors(1) == 1, "tau 1");
    ASSERT_TRUE(&ctx, ml_sigma_sum(36) == 91, "sigma 36");
    ASSERT_TRUE(&ctx, ml_sigma_sum(6) == 12, "sigma 6 perfect");
    ASSERT_TRUE(&ctx, ml_frobenius2(4, 7) == 17, "Frob 4,7");
    ASSERT_TRUE(&ctx, ml_frobenius2(4, 6) == 0, "Frob non-coprime 0");
    ASSERT_TRUE(&ctx, ml_frobenius_count2(4, 7) == 9, "Frob count 4,7");

    /* Tier 3: series (vs Python/mpmath truth). */
    ASSERT_NEAR(&ctx, ml_harmonic_number(10), 2.9289682539682538, 1e-12, "H10");
    ASSERT_NEAR(&ctx, ml_harmonic_number(1), 1.0, 0.0, "H1");
    ASSERT_TRUE(&ctx, ml_isnan(ml_harmonic_number(0)), "H0 NaN");
    ASSERT_NEAR(&ctx, ml_basel_partial(1000), 1.6439345666815599, 1e-12, "Basel1000");
    ASSERT_NEAR(&ctx, ml_faulhaber_sum(100, 4), 2050333330.0, 1e-3, "Faulhaber 100^4");
    ASSERT_NEAR(&ctx, ml_faulhaber_sum(10, 1), 55.0, 1e-9, "Faulhaber 10^1");

    /* Tier 3: Wythoff / Grundy / inequalities / Vieta. */
    {
        uint64_t a = 0, b = 0;
        ml_status_t st = ml_wythoff_cold(5, &a, &b);
        ASSERT_TRUE(&ctx, st == ML_SUCCESS && a == 8 && b == 13, "Wythoff k=5");
        st = ml_wythoff_cold(10, &a, &b);
        ASSERT_TRUE(&ctx, st == ML_SUCCESS && a == 16 && b == 26, "Wythoff k=10");
        ASSERT_TRUE(&ctx, ml_wythoff_cold(0, &a, &b) == ML_ERR_INVALID_ARG, "Wythoff 0 invalid");
    }
    {
        int mv[] = {1, 2, 3};
        ASSERT_TRUE(&ctx, ml_subtract_game_grundy(0, mv, 3) == 0, "Grundy 0");
        ASSERT_TRUE(&ctx, ml_subtract_game_grundy(4, mv, 3) == 0, "Grundy 4 P-pos");
        ASSERT_TRUE(&ctx, ml_subtract_game_grundy(5, mv, 3) == 1, "Grundy 5");
        ASSERT_TRUE(&ctx, ml_subtract_game_grundy(5000, mv, 3) == -1, "Grundy cap");
    }
    ASSERT_TRUE(&ctx, ml_schur_deg1(1, 2, 3) == 1, "Schur holds");
    ASSERT_TRUE(&ctx, ml_schur_deg1(0, 0, 0) == 1, "Schur zeros");
    ASSERT_TRUE(&ctx, ml_schur_deg1(-1, 2, 3) == 0, "Schur negative invalid");
    {
        double d[] = {1, 2, 3, 4};
        ASSERT_TRUE(&ctx, ml_amgm_check(d, 4) == 1, "AMGM holds");
        ASSERT_TRUE(&ctx, ml_amgm_check(d, 0) == 0, "AMGM invalid");
    }
    {
        double s1 = 0, s2 = 0, s3 = 0;
        /* x^3-6x^2+11x-6 = (x-1)(x-2)(x-3). */
        ml_status_t st = ml_vieta_cubic_sums(1, -6, 11, -6, &s1, &s2, &s3);
        ASSERT_TRUE(&ctx, st == ML_SUCCESS, "Vieta status");
        ASSERT_NEAR(&ctx, s1, 6.0, 1e-12, "Vieta s1");
        ASSERT_NEAR(&ctx, s2, 11.0, 1e-12, "Vieta s2");
        ASSERT_NEAR(&ctx, s3, 6.0, 1e-12, "Vieta s3");
        ASSERT_TRUE(&ctx, ml_vieta_cubic_sums(0, 1, 2, 3, &s1, &s2, &s3) == ML_ERR_INVALID_ARG, "Vieta a=0");
    }

    return ml_test_summary(&ctx);
}
