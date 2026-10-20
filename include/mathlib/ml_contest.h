#ifndef MATHLIB_ML_CONTEST_H
#define MATHLIB_ML_CONTEST_H

#include <stdint.h>
#include "ml_compiler.h"
#include "ml_core.h"
#include "ml_types.h"
#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Contest mathematics: recreational -> olympiad -> Putnam.
 *
 * Extends the existing contest core (ml_combinatorics: factorial/nPr/nCr/
 * catalan/fib/derange; ml_numbertheory: legendre/jacobi/mobius/crt/prime_pi;
 * ml_numerical: arith/geom sums, sum_k/k2/k3, nim, majorizes; ml_analytic_nt:
 * partition/theta) with the classic missing staples, tiered low -> high.
 *
 * Exactness contract (matches house style):
 * - u64 functions are EXACT (0 ULP); UINT64_MAX is the overflow/invalid
 *   sentinel (hanoi(64) would equal it, so it is capped, never returned
 *   as a value).
 * - int verdicts are 1/0; invalid input fails closed to 0 (nim_win,
 *   majorizes convention), documented per function.
 * - double sums are exactness-guarded like ml_sum_k: NaN beyond the range
 *   where the closed form round-trips through binary64.
 * - All scratch is stack-local or heap-per-call (no shared mutable state);
 *   one workspace/thread rule from ml_tensor.h applies where workspaces
 *   are used (not needed here: caps keep every frame < 20KB).
 * - No oracle-tier ULP certification claimed (PRECISION_CONTRACT rule 3:
 *   edge-tier exact + regression-guard tier only).
 * ========================================================================== */

/* ---- Tier 0: recreational (AMC 8/10 kit) ---- */
ML_API uint64_t ml_digit_sum_u64 (uint64_t n);
ML_API uint64_t ml_digital_root (uint64_t n);
ML_API int ml_is_palindrome_u64 (uint64_t n);
ML_API uint64_t ml_collatz_steps (uint64_t n);
ML_API int ml_magic_square_check (const long long *a, int n);

/* ---- Tier 1: school / AMC (figurate, parity, recurrences, games) ---- */
ML_API uint64_t ml_triangular (uint64_t n);
ML_API uint64_t ml_pentagonal (uint64_t n);
ML_API uint64_t ml_hexagonal (uint64_t n);
ML_API int ml_is_triangular (uint64_t m);
ML_API uint64_t ml_double_factorial (long long n);
ML_API int ml_binom_odd (uint64_t n, uint64_t k);
ML_API uint64_t ml_stars_bars (uint64_t stars, uint64_t bars);
ML_API uint64_t ml_hanoi_moves (int n);
ML_API uint64_t ml_josephus (uint64_t n, uint64_t k);
ML_API uint64_t ml_lucas_number (int n);
ML_API uint64_t ml_pell_number (int n);
ML_API uint64_t ml_tribonacci (int n);

/* ---- Tier 2: AIME / olympiad (partitions, primes, congruences) ---- */
ML_API uint64_t ml_bell_number (int n);
ML_API uint64_t ml_stirling2 (int n, int k);
ML_API uint64_t ml_stirling1 (int n, int k);
ML_API uint64_t ml_eulerian (int n, int m);
ML_API uint64_t ml_lah_number (int n, int k);
ML_API uint64_t ml_motzkin (int n);
ML_API uint64_t ml_schroder_large (int n);
ML_API uint64_t ml_narayana (int n, int k);
ML_API uint64_t ml_part_distinct (int n);
ML_API uint64_t ml_goldbach_partitions (uint64_t n);
ML_API int ml_is_sophie_germain (uint64_t p);
ML_API int ml_lucas_lehmer (int p);
ML_API uint64_t ml_modsqrt_tonelli (uint64_t a, uint64_t p);
ML_API uint64_t ml_tau_divisors (uint64_t n);
ML_API uint64_t ml_sigma_sum (uint64_t n);
ML_API uint64_t ml_frobenius2 (uint64_t a, uint64_t b);
ML_API uint64_t ml_frobenius_count2 (uint64_t a, uint64_t b);

/* ---- Tier 3: Putnam / college (series, games, inequalities, vieta) ---- */
ML_API double ml_harmonic_number (long long n);
ML_API double ml_basel_partial (long long n);
ML_API double ml_faulhaber_sum (long long n, int p);
ML_API ml_status_t ml_wythoff_cold (uint64_t k, uint64_t *a, uint64_t *b);
ML_API int ml_subtract_game_grundy (uint64_t n, const int *moves, int nmoves);
ML_API int ml_schur_deg1 (double a, double b, double c);
ML_API int ml_amgm_check (const double *data, int n);
ML_API ml_status_t ml_vieta_cubic_sums (double a, double b, double c, double d,
                                        double *s1, double *s2, double *s3);

#ifdef __cplusplus
}
#endif
#endif /* MATHLIB_ML_CONTEST_H */
