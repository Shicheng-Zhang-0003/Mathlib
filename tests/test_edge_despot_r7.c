/* Audit regressions: libm parity, reflection tails, and fixed-capacity bigint. */
#include "test_harness.h"
#include "ml_bigint.h"
#include "ml_graph.h"
#include "ml_integral.h"
#include "ml_math_extra.h"
#include "ml_special.h"
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
static uint64_t audit_next (uint64_t *state) {
    uint64_t x = *state;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    *state = x;
    return x * 2685821657736338717ULL;
}
int main (void) {
    ml_test_ctx_t ctx;
    ml_test_init (&ctx, "Despot R7 audit regressions");

    /* remquo must return quotient bits without narrowing the full quotient. */
    {
        int q = 0;
        double r = ml_remquo (0x1p100, 3.0, &q);
        ASSERT_NEAR (&ctx, r, 1.0, 0.0, "large remquo remainder");
        ASSERT_TRUE (&ctx, q == 85, "large remquo low seven bits");
        r = ml_remquo (-0x1p100, 3.0, &q);
        ASSERT_NEAR (&ctx, r, -1.0, 0.0, "negative large remquo remainder");
        ASSERT_TRUE (&ctx, q == -85, "negative remquo quotient sign");
        r = ml_remquo (7.0, 2.0, &q);
        ASSERT_NEAR (&ctx, r, -1.0, 0.0, "remquo tie remainder");
        ASSERT_TRUE (&ctx, q == 4, "remquo nearest-even quotient");
    }

    /* Compare low quotient bits and exact remainders against the platform
     * remquo over exponent-spread finite binary64 inputs. */
    {
        uint64_t state = 0xD1B54A32D192ED03ULL;
        int i;
        for (i = 0; i < 128; i++) {
            uint64_t rx = audit_next (&state), ry = audit_next (&state);
            int ex = (int) (rx % 2097u) - 1074;
            int ey = (int) (ry % 2097u) - 1074;
            double mx = 1.0 + (double) ((rx >> 12) & 0xFFFFFULL) / 1048576.0;
            double my = 1.0 + (double) ((ry >> 12) & 0xFFFFFULL) / 1048576.0;
            double x = ldexp (mx, ex), y = ldexp (my, ey);
            int qlib = 0, qsys = 0;
            double rlib, rsys;
            if (rx >> 63) x = -x;
            if (ry >> 63) y = -y;
            rlib = ml_remquo (x, y, &qlib);
            rsys = remquo (x, y, &qsys);
            ASSERT_TRUE (&ctx, rlib == rsys && ml_signbit (rlib) == ml_signbit (rsys),
                         "remquo randomized remainder matches libm");
            ASSERT_TRUE (&ctx, qlib % 8 == qsys % 8,
                         "remquo randomized low quotient bits match libm");
        }
    }

    /* C99 modf propagates NaN to both results, and integer casts stay defined. */
    {
        double ip = 0.0;
        double f = ml_modf (NAN, &ip);
        ASSERT_TRUE (&ctx, ml_isnan (f) && ml_isnan (ip), "modf NaN propagation");
        f = ml_modf (-3.0, &ip);
        ASSERT_TRUE (&ctx, f == 0.0 && ml_signbit (f) && ml_signbit (ip),
                     "modf negative integer signed fractional zero");
        f = ml_modf (-0.0, &ip);
        ASSERT_TRUE (&ctx, f == 0.0 && ml_signbit (f) && ml_signbit (ip),
                     "modf preserves negative zero");
        ASSERT_TRUE (&ctx, ml_llrint (DBL_MAX) == LLONG_MIN, "llrint overflow sentinel");
        ASSERT_TRUE (&ctx, ml_llround (DBL_MAX) == LLONG_MIN, "llround overflow sentinel");
        ASSERT_TRUE (&ctx, ml_lrint (0.5) == 0L && ml_lround (0.5) == 1L,
                     "rint and round modes");
    }

    /* sign(Gamma(x)) is the sign of sin(pi*x) on negative non-poles. */
    {
        int sign = 0;
        (void) ml_lgamma_r (-0.5, &sign);
        ASSERT_TRUE (&ctx, sign == -1, "lgamma_r sign at -1/2");
        (void) ml_lgamma_r (-1.5, &sign);
        ASSERT_TRUE (&ctx, sign == 1, "lgamma_r sign at -3/2");
        (void) ml_lgamma_r (-2.0, &sign);
        ASSERT_TRUE (&ctx, sign == 0, "lgamma_r pole sign is undefined");
    }

    /* Large negative binary64 inputs are even integers and trivial zeta zeros. */
    ASSERT_TRUE (&ctx, ml_zeta (-0x1p53) == 0.0, "zeta zero at -2^53");
    ASSERT_TRUE (&ctx, ml_zeta (-DBL_MAX) == 0.0, "zeta zero at -DBL_MAX");

    /* Empty graph inputs do not rely on qsort accepting a null base pointer. */
    {
        int dummy_edge = 0;
        double dummy_weight = 0.0, total = -1.0;
        ASSERT_TRUE (&ctx, ml_kruskal_mst (1, &dummy_edge, &dummy_edge, &dummy_weight, 0, &total) == ML_SUCCESS && total == 0.0,
                     "empty one-vertex MST");
    }

    /* Newly added special functions: identities and finite boundary limits. */
    {
        double k = ml_ellip_k (ml_sqrt (0.5));
        double u = 2.0 * k - 0.1;
        double sn = ml_jacobi_sn (u, 0.5);
        double cn = ml_jacobi_cn (u, 0.5);
        double dn = ml_jacobi_dn (u, 0.5);
        double sn_near_one = ml_jacobi_sn (10.0, 0.99);
        double cn_near_one = ml_jacobi_cn (10.0, 0.99);
        double dn_near_one = ml_jacobi_dn (10.0, 0.99);
        ASSERT_NEAR (&ctx, ml_bessel_jv (-1.0, 1.0), -0.4400505857449335, 2e-13,
                     "negative integer order J reflection");
        ASSERT_NEAR (&ctx, ml_bessel_iv (-1.0, 1.0), 0.5651591039924850, 2e-13,
                     "negative integer order I reflection");
        ASSERT_NEAR (&ctx, ml_bessel_jv (-1.5, 1.0), -1.1024955751601792, 2e-13,
                     "negative noninteger J gamma sign");
        ASSERT_NEAR (&ctx, ml_bessel_iv (-1.5, 1.0), -0.2935253263474798, 2e-13,
                     "negative noninteger I gamma sign");
        ASSERT_NEAR (&ctx, ml_bessel_jv (0.1, 17.9), -0.06072988900411563, 3e-12,
                     "general-order J near series cancellation limit");
        ASSERT_NEAR (&ctx, ml_sph_bessel_k (0, 1.0), 0.5778636748954609, 2e-13,
                     "spherical k0 normalization");
        ASSERT_NEAR (&ctx, ml_bessel_jv (20.0, 18.0), 0.0673059474374060, 2e-14,
                     "integer-order Miller recurrence near turning point");
        ASSERT_NEAR (&ctx, ml_bessel_yv (20.0, 18.0), -0.5542283743188141, 2e-13,
                     "high integer-order Y recurrence");
        ASSERT_NEAR (&ctx, ml_bessel_iv (20.0, 18.0), 178.8238343018977, 2e-11,
                     "positive-term high-order I series");
        ASSERT_TRUE (&ctx, ml_fabs (ml_bessel_iv (20.0, 700.0) / 1.1492421734080173e302 - 1.0) < 5e-13,
                     "large finite I value uses converged positive series");
        ASSERT_NEAR (&ctx, ml_bessel_kv (20.0, 18.0), 1.039001642536835e-4, 2e-16,
                     "high integer-order K recurrence");
        ASSERT_NEAR (&ctx, ml_bessel_kv (0.1, 17.9), 4.953677122532237e-9, 2e-20,
                     "noninteger K integral avoids cancellation");
        ASSERT_NEAR (&ctx, ml_bessel_kv (20.1, 30.0), 1.309020411562299e-11, 2e-23,
                     "noninteger high-order K integral");
        ASSERT_TRUE (&ctx, ml_isnan (ml_bessel_yv (1.000000000001, 1.0)) && ml_isnan (ml_bessel_kv (1.000000000001, 1.0)),
                     "near-integer Y/K cancellation rejected explicitly");
        ASSERT_TRUE (&ctx, ml_isnan (ml_bessel_jv (20.1, 30.0)) && ml_isnan (ml_bessel_yv (20.1, 30.0)),
                     "unsupported fixed-order transition rejected explicitly");
        ASSERT_TRUE (&ctx, ml_fabs (ml_polygamma (100, 12.0) / (-9.393593619882782e48) - 1.0) < 2e-13,
                     "high-order polygamma uses convergent tail sum");
        ASSERT_NEAR (&ctx, ml_airy_bi (0.0), 0.6149266274460007, 2e-15,
                     "Airy Bi normalization at zero");
        ASSERT_NEAR (&ctx, ml_airy_bi (-10.0), -0.3146798296438386, 2e-13,
                     "Airy Bi oscillatory branch");
        ASSERT_NEAR (&ctx, ml_airy_ai_prime (10.0), -3.5206336767389237e-10, 2e-22,
                     "Airy Ai derivative decaying branch");
        ASSERT_NEAR (&ctx, ml_airy_ai_prime (-10.0), 0.9962650441327900, 5e-13,
                     "Airy Ai derivative oscillatory branch");
        ASSERT_NEAR (&ctx, ml_airy_bi_prime (-10.0), 0.1194141133999092, 1e-12,
                     "Airy Bi derivative oscillatory branch");
        ASSERT_NEAR (&ctx, ml_dawson (3.5), 0.1496215930807565, 2e-13,
                     "Dawson series and asymptotic seam");
        ASSERT_TRUE (&ctx, ml_si (ml_make_inf (0)) == 1.5707963267948966,
                     "Si positive-infinity limit");
        ASSERT_TRUE (&ctx, ml_ci (ml_make_inf (0)) == 0.0,
                     "Ci positive-infinity limit");
        ASSERT_NEAR (&ctx, ml_fresnel_c (1.0), 0.7798934003768230, 2e-13,
                     "Fresnel C at one");
        ASSERT_NEAR (&ctx, ml_fresnel_s (1.0), 0.4382591473903548, 2e-13,
                     "Fresnel S at one");
        ASSERT_NEAR (&ctx, ml_fresnel_c (4.0), 0.4984260330381776, 2e-10,
                     "Fresnel C asymptotic branch");
        ASSERT_NEAR (&ctx, ml_fresnel_s (4.0), 0.4205157542469284, 2e-10,
                     "Fresnel S asymptotic branch");
        ASSERT_TRUE (&ctx, ml_fresnel_c (ml_make_inf (0)) == 0.5 && ml_fresnel_s (ml_make_inf (1)) == -0.5,
                     "Fresnel signed infinite limits");
        ASSERT_TRUE (&ctx, ml_sph_bessel_j (0, -1.0) == ml_sph_bessel_j (0, 1.0) && ml_sph_bessel_j (1, -1.0) == -ml_sph_bessel_j (1, 1.0) && ml_sph_bessel_i (1, -1.0) == -ml_sph_bessel_i (1, 1.0),
                     "spherical j/i parity on negative arguments");
        ASSERT_TRUE (&ctx, ml_sph_bessel_y (0, -1.0) == -ml_sph_bessel_y (0, 1.0) && ml_sph_bessel_y (1, -1.0) == ml_sph_bessel_y (1, 1.0),
                     "spherical y parity on negative arguments");
        ASSERT_TRUE (&ctx, ml_isinf (ml_sph_bessel_i (0, ml_make_inf (0))) && ml_isinf (ml_sph_bessel_i (1, ml_make_inf (1))) && ml_signbit (ml_sph_bessel_i (1, ml_make_inf (1))),
                     "spherical i infinite limits preserve parity");
        ASSERT_TRUE (&ctx, ml_isnan (ml_ci (-1.0)) && ml_isnan (ml_chi (-1.0)),
                     "Ci and Chi reject negative real branch inputs");
        ASSERT_TRUE (&ctx, cn < 0.0 && ml_fabs (sn * sn + cn * cn - 1.0) < 2e-14 && ml_fabs (dn * dn + 0.5 * sn * sn - 1.0) < 2e-14,
                     "Jacobi elliptic quadrant and identities");
        ASSERT_NEAR (&ctx, sn_near_one, -0.9914207448606821, 2e-14,
                     "Jacobi bisection at large u near m=1: sn");
        ASSERT_NEAR (&ctx, cn_near_one, -0.1307092447376634, 2e-14,
                     "Jacobi bisection at large u near m=1: cn");
        ASSERT_NEAR (&ctx, dn_near_one, 0.1640550444006264, 2e-14,
                     "Jacobi bisection at large u near m=1: dn");
        ASSERT_NEAR (&ctx, ml_struve_h0 (50.0), -0.0853376748261190, 2e-14,
                     "Struve H0 oscillatory large-argument integral");
        ASSERT_NEAR (&ctx, ml_struve_h0 (100.0), -0.07087875168964734, 2e-14,
                     "Struve H0 large-argument accuracy");
        ASSERT_TRUE (&ctx, ml_isnan (ml_struve_h0 (10001.0)),
                     "Struve H0 documented work bound");
        ASSERT_TRUE (&ctx, ml_isnan (ml_bessel_yv (5000.0, 1.0)) && ml_isnan (ml_assoc_legendre (4, INT_MIN, 0.5)),
                     "unsupported order and INT_MIN rejected");
    }

    /* The DD gamma path must be used consistently above x=8. */
    ASSERT_TRUE (&ctx, ml_gamma_new (150.125) == 7.122871194429566e260,
                 "gamma correctly rounded at 150.125");

    /* Small-composite certainty, large probable-prime rounds, and RNG bounds. */
    {
        ml_bigint_t n, r;
        uint64_t state = 0u;
        ml_bigint_init (&n);
        (void) ml_bigint_set_u64 (&n, 41u * 41u);
        ASSERT_TRUE (&ctx, !ml_bigint_is_probable_prime (&n, 8, &state),
                     "41 squared is composite");
        (void) ml_bigint_set_u64 (&n, 41u);
        ASSERT_TRUE (&ctx, ml_bigint_is_probable_prime (&n, 8, &state),
                     "41 is prime");
        ASSERT_TRUE (&ctx, ml_bigint_rand (&r, UINT_MAX, &state) == ML_ERR_INVALID_ARG,
                     "random bit count overflow rejected");
        ASSERT_TRUE (&ctx, ml_bigint_rand (&r, 8u, &state) == ML_SUCCESS && ml_bigint_bitlen (&r) <= 8 && state != 0u,
                     "zero seed has bounded nonzero RNG stream");
    }

    /* Text conversion accepts C whitespace at the edges and never commits a
     * partial parse when malformed input is encountered. */
    {
        ml_bigint_t n;
        (void) ml_bigint_set_i64 (&n, 17);
        ASSERT_TRUE (&ctx, ml_bigint_from_str (&n, "\r\n -42\t", 10) == ML_SUCCESS && ml_bigint_to_double (&n) == -42.0,
                     "parser leading and trailing whitespace");
        ASSERT_TRUE (&ctx, ml_bigint_from_str (&n, "12 \n 3", 10) == ML_ERR_INVALID_ARG && ml_bigint_to_double (&n) == -42.0,
                     "parser rejects digits after whitespace transactionally");
        ASSERT_TRUE (&ctx, ml_bigint_from_str (&n, "  +", 10) == ML_ERR_INVALID_ARG,
                     "parser rejects missing digits");
    }

    /* Small-value differential properties cover signed arithmetic and the
     * Euclidean-modulo convention against native integer arithmetic. */
    {
        uint64_t state = 0xA0761D6478BD642FULL;
        int i;
        for (i = 0; i < 128; i++) {
            int64_t av = (int64_t) (audit_next (&state) % 2000001u) - 1000000;
            int64_t bv = (int64_t) (audit_next (&state) % 2000001u) - 1000000;
            ml_bigint_t a, b, r;
            int64_t expected_mod;
            if (bv == 0) bv = 1;
            (void) ml_bigint_set_i64 (&a, av);
            (void) ml_bigint_set_i64 (&b, bv);
            ASSERT_TRUE (&ctx, ml_bigint_add (&r, &a, &b) == ML_SUCCESS && ml_bigint_to_double (&r) == (double) (av + bv),
                         "bigint randomized addition");
            ASSERT_TRUE (&ctx, ml_bigint_sub (&r, &a, &b) == ML_SUCCESS && ml_bigint_to_double (&r) == (double) (av - bv),
                         "bigint randomized subtraction");
            ASSERT_TRUE (&ctx, ml_bigint_mul (&r, &a, &b) == ML_SUCCESS && ml_bigint_to_double (&r) == (double) (av * bv),
                         "bigint randomized multiplication");
            ASSERT_TRUE (&ctx, ml_bigint_div (&r, NULL, &a, &b) == ML_SUCCESS && ml_bigint_to_double (&r) == (double) (av / bv),
                         "bigint randomized truncating division");
            expected_mod = av % bv;
            if (expected_mod < 0) expected_mod += (bv < 0) ? -bv : bv;
            ASSERT_TRUE (&ctx, ml_bigint_mod (&r, &a, &b) == ML_SUCCESS && ml_bigint_to_double (&r) == (double) expected_mod,
                         "bigint randomized Euclidean modulo");
        }
    }

    /* A full-capacity value can be shifted by zero and divided without a
     * spurious overflow from a temporary carry limb. */
    {
        ml_bigint_t a, shifted, one, q, rem, keep;
        ml_bigint_init (&a);
        a.limb [ML_BIGINT_LIMBS - 1] = 0x80000000u;
        a.n = ML_BIGINT_LIMBS;
        ASSERT_TRUE (&ctx, ml_bigint_shl (&shifted, &a, 0u) == ML_SUCCESS,
                     "full-width zero shift succeeds");
        ASSERT_TRUE (&ctx, ml_bigint_cmp (&shifted, &a) == 0,
                     "full-width zero shift preserves value");
        (void) ml_bigint_set_u64 (&one, 1u);
        ASSERT_TRUE (&ctx, ml_bigint_div (&q, &rem, &a, &one) == ML_SUCCESS,
                     "full-width division by one");
        ASSERT_TRUE (&ctx, ml_bigint_cmp (&q, &a) == 0 && ml_bigint_is_zero (&rem),
                     "full-width division identity");
        (void) ml_bigint_set_u64 (&keep, 7u);
        ASSERT_TRUE (&ctx, ml_bigint_shl (&keep, &a, 1u) == ML_ERR_OVERFLOW,
                     "full-width overflowing shift rejected");
        ASSERT_TRUE (&ctx, ml_bigint_to_double (&keep) == 7.0,
                     "overflow leaves output unchanged");
        ASSERT_TRUE (&ctx, ml_bigint_add (&shifted, &a, &one) == ML_SUCCESS,
                     "full-width addition that fits succeeds");
        keep.n = 1;
        keep.limb [0] = 0u;
        ASSERT_TRUE (&ctx, ml_bigint_add (&shifted, &a, &keep) == ML_ERR_INVALID_ARG,
                     "noncanonical operand rejected");
    }

    /* Euclidean modulo is independent of divisor sign; modular power handles
     * identity exponent and products wider than the result capacity. */
    {
        ml_bigint_t a, b, r, exp, mod, base, expected;
        (void) ml_bigint_set_i64 (&a, -7);
        (void) ml_bigint_set_i64 (&b, -3);
        ASSERT_TRUE (&ctx, ml_bigint_mod (&r, &a, &b) == ML_SUCCESS && ml_bigint_to_double (&r) == 2.0,
                     "negative divisor Euclidean modulo");
        (void) ml_bigint_set_u64 (&a, 19u);
        (void) ml_bigint_set_u64 (&exp, 0u);
        (void) ml_bigint_set_u64 (&mod, 1u);
        ASSERT_TRUE (&ctx, ml_bigint_powmod (&r, &a, &exp, &mod) == ML_SUCCESS && ml_bigint_is_zero (&r), "zero exponent modulo one");

        ml_bigint_init (&mod);
        for (int i = 0; i < ML_BIGINT_LIMBS; i++) mod.limb [i] = UINT32_MAX;
        mod.n = ML_BIGINT_LIMBS;
        ml_bigint_init (&base);
        base.limb [ML_BIGINT_LIMBS / 2] = 1u;
        base.n = ML_BIGINT_LIMBS / 2 + 1;
        (void) ml_bigint_set_u64 (&exp, 2u);
        (void) ml_bigint_set_u64 (&expected, 1u);
        ASSERT_TRUE (&ctx, ml_bigint_powmod (&r, &base, &exp, &mod) == ML_SUCCESS && ml_bigint_cmp (&r, &expected) == 0,
                     "full-width modular square avoids wide product overflow");
    }

    /* Signed modular inverse and divide-before-multiply LCM preserve valid
     * values even when the direct product would overflow the fixed capacity. */
    {
        ml_bigint_t a, m, r, expected, x, y;
        (void) ml_bigint_set_i64 (&a, -3);
        (void) ml_bigint_set_u64 (&m, 11u);
        (void) ml_bigint_set_u64 (&expected, 7u);
        ASSERT_TRUE (&ctx, ml_bigint_modinv (&r, &a, &m) == ML_SUCCESS && ml_bigint_cmp (&r, &expected) == 0,
                     "negative modular inverse");
        ml_bigint_init (&x);
        x.limb [ML_BIGINT_LIMBS - 1] = 0x80000000u;
        x.n = ML_BIGINT_LIMBS;
        y = x;
        ASSERT_TRUE (&ctx, ml_bigint_lcm (&r, &x, &y) == ML_SUCCESS && ml_bigint_cmp (&r, &x) == 0,
                     "LCM divides before a product that would overflow");
    }

    return ml_test_summary (&ctx);
}
