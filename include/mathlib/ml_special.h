#ifndef MATHLIB_ML_SPECIAL_H
#define MATHLIB_ML_SPECIAL_H

#include "ml_compiler.h"
#include "ml_core.h"
#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * v12R3: extended special functions.
 *
 * Bessel (general real order), spherical Bessel, Airy Bi + derivatives,
 * exponential/logarithmic integrals, Fresnel, Dawson, scaled complementary
 * error, polygamma, orthogonal polynomials (associated Legendre, Gegenbauer,
 * Jacobi, Zernike), real spherical harmonics, polylogarithm, Lerch
 * transcendent, complete elliptic integral of the third kind, Jacobi
 * elliptic functions, Struve. This module is experimental and has no
 * oracle-tier accuracy certification. Integer orders/degrees are bounded to
 * 4096; polygamma orders are bounded to 100 (n>0 currently requires x>0).
 * Series-based polylog/Lerch calls return NaN when their fixed iteration
 * budget does not meet the convergence test. For noninteger Bessel orders,
 * fixed-order asymptotics are used only for x>=max(18,4*nu^2); unsupported
 * high-order transition cases return NaN. Yv/Kv return NaN when a noninteger
 * order is within 1e-7 of an integer, where their connection formulas lose
 * precision; exact integer orders use recurrence. ellip_pi supports
 * 0<=k<1,n<1. Struve H0 uses an oscillatory integral on |x|<=10000 and
 * returns NaN outside that finite range.
 * ========================================================================== */

/* Bessel functions of general real order (x >= 0 for Jv/Iv; x > 0 for Yv/Kv) */
ML_API double ml_bessel_jv (double nu, double x);
ML_API double ml_bessel_yv (double nu, double x);
ML_API double ml_bessel_iv (double nu, double x);
ML_API double ml_bessel_kv (double nu, double x);
ML_API double ml_bessel_jn (int n, double x);
ML_API double ml_bessel_yn (int n, double x);

/* Spherical j/i/y obey their integer-degree parity on negative x; spherical
 * k is defined here for x>=0 (with its singular limit at zero). */
ML_API double ml_sph_bessel_j (int l, double x);
ML_API double ml_sph_bessel_y (int l, double x);
ML_API double ml_sph_bessel_i (int l, double x);
ML_API double ml_sph_bessel_k (int l, double x);

/* Airy functions */
ML_API double ml_airy_bi (double x);
ML_API double ml_airy_ai_prime (double x);
ML_API double ml_airy_bi_prime (double x);

/* E1 is real-valued for x>=0; Ci and Chi require x>=0; Si and Shi accept
 * either sign of x. */
ML_API double ml_expint_e1 (double x); /* E1(x) */
ML_API double ml_si (double x);        /* Si(x) */
ML_API double ml_ci (double x);        /* Ci(x) */
ML_API double ml_shi (double x);       /* Shi(x) */
ML_API double ml_chi (double x);       /* Chi(x) */

/* Fresnel, Dawson, scaled erfc */
ML_API double ml_fresnel_s (double x);
ML_API double ml_fresnel_c (double x);
ML_API double ml_dawson (double x);
ML_API double ml_erfcx (double x);

/* Polygamma: n in [0,100]; n>0 currently requires x>0. */
ML_API double ml_polygamma (int n, double x);
ML_API double ml_trigamma (double x);

/* Orthogonal polynomials */
ML_API double ml_assoc_legendre (int l, int m, double x);
ML_API double ml_spherical_harmonic (int l, int m, double theta, double phi);
ML_API double ml_gegenbauer (int n, double alpha, double x);
ML_API double ml_jacobi_p (int n, double a, double b, double x);
ML_API double ml_zernike_r (int n, int m, double rho);

/* Polylogarithm and Lerch transcendent */
ML_API double ml_polylog (int s, double x);
ML_API double ml_lerch_phi (double z, double s, double a);

/* Complete elliptic integral of the third kind (0<=k<1,n<1), and Jacobi
 * elliptic functions (parameter m = k^2). */
ML_API double ml_ellip_pi (double n, double k);
ML_API double ml_jacobi_sn (double u, double m);
ML_API double ml_jacobi_cn (double u, double m);
ML_API double ml_jacobi_dn (double u, double m);

/* Struve */
ML_API double ml_struve_h0 (double x);

#ifdef __cplusplus
}
#endif
#endif /* MATHLIB_ML_SPECIAL_H */
