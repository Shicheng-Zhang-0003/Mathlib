#include "test_harness.h"
#include "ml_polyq.h"

#include <string.h>
static int mpq_equals (const ml_mpq_t *q, const char *num, const char *den) {
    char nb [128], db [128];
    return ml_mpz_to_string (&q->numerator, nb, sizeof (nb), 10u) >= 0 &&
           ml_mpz_to_string (&q->denominator, db, sizeof (db), 10u) >= 0 &&
           strcmp (nb, num) == 0 && strcmp (db, den) == 0;
}
static ml_status_t set_coeff_i64 (ml_polyq_t *p, size_t degree, int64_t value) {
    ml_mpq_t q;
    ml_status_t st;
    ml_mpq_init (&q);
    st = ml_mpq_set_i64 (&q, value);
    if (st == ML_SUCCESS) st = ml_polyq_set_coefficient (p, degree, &q);
    ml_mpq_clear (&q);
    return st;
}
int main (void) {
    ml_test_ctx_t ctx;
    ml_test_init (&ctx, "Exact symbolic polynomial foundation");

    {
        ml_polyq_t x, f, g, quotient, remainder, product, derivative, gcd;
        ml_mpq_t at_two, value;
        int certificate = 0, equal = 0;
        ml_polyq_init (&x);
        ml_polyq_init (&f);
        ml_polyq_init (&g);
        ml_polyq_init (&quotient);
        ml_polyq_init (&remainder);
        ml_polyq_init (&product);
        ml_polyq_init (&derivative);
        ml_polyq_init (&gcd);
        ml_mpq_init (&at_two);
        ml_mpq_init (&value);

        ASSERT_TRUE (&ctx, set_coeff_i64 (&x, 1u, 1) == ML_SUCCESS && ml_polyq_degree (&x) == 1,
                     "construct the exact indeterminate x");
        ASSERT_TRUE (&ctx, set_coeff_i64 (&f, 0u, -1) == ML_SUCCESS && set_coeff_i64 (&f, 2u, 1) == ML_SUCCESS,
                     "construct x squared minus one");
        ASSERT_TRUE (&ctx, set_coeff_i64 (&g, 0u, -1) == ML_SUCCESS && set_coeff_i64 (&g, 1u, 1) == ML_SUCCESS,
                     "construct x minus one");
        ASSERT_TRUE (&ctx, ml_polyq_mul (&product, &g, &x) == ML_SUCCESS && ml_polyq_add (&product, &product, &g) == ML_SUCCESS && ml_polyq_equal (&product, &f, &equal) == ML_SUCCESS && equal,
                     "exact polynomial add and multiply recover x squared minus one");

        ASSERT_TRUE (&ctx, ml_polyq_divrem (&quotient, &remainder, &f, &g) == ML_SUCCESS && ml_polyq_degree (&quotient) == 1 && ml_polyq_is_zero (&remainder),
                     "exact polynomial long division");
        ASSERT_TRUE (&ctx, ml_polyq_verify_division (&f, &g, &quotient, &remainder, &certificate) == ML_SUCCESS && certificate,
                     "independent exact division-certificate verification");
        ASSERT_TRUE (&ctx, ml_polyq_get_coefficient (&quotient, 0u, &value) == ML_SUCCESS && mpq_equals (&value, "1", "1") && ml_polyq_get_coefficient (&quotient, 1u, &value) == ML_SUCCESS && mpq_equals (&value, "1", "1"),
                     "quotient coefficients are exact and reduced");

        ASSERT_TRUE (&ctx, set_coeff_i64 (&remainder, 0u, 1) == ML_SUCCESS && ml_polyq_verify_division (&f, &g, &quotient, &remainder, &certificate) == ML_SUCCESS && !certificate,
                     "certificate rejects a forged remainder");

        ASSERT_TRUE (&ctx, ml_polyq_derivative (&derivative, &f) == ML_SUCCESS && ml_polyq_degree (&derivative) == 1 && ml_polyq_get_coefficient (&derivative, 1u, &value) == ML_SUCCESS && mpq_equals (&value, "2", "1"),
                     "symbolic derivative is exact");
        (void) ml_mpq_set_i64 (&at_two, 2);
        ASSERT_TRUE (&ctx, ml_polyq_evaluate (&value, &f, &at_two) == ML_SUCCESS && mpq_equals (&value, "3", "1"),
                     "Horner evaluation remains exact over rationals");
        ASSERT_TRUE (&ctx, ml_polyq_gcd (&gcd, &f, &g) == ML_SUCCESS && ml_polyq_equal (&gcd, &g, &equal) == ML_SUCCESS && equal,
                     "Euclidean polynomial gcd returns monic common factor");
        ASSERT_TRUE (&ctx, ml_polyq_divrem (&quotient, &remainder, &f, &x) == ML_SUCCESS && ml_polyq_verify_division (&f, &x, &quotient, &remainder, &certificate) == ML_SUCCESS && certificate,
                     "division certificate covers nonzero remainder");
        (void) set_coeff_i64 (&remainder, 0u, 0);
        ASSERT_TRUE (&ctx, ml_polyq_divrem (&quotient, &remainder, &f, &remainder) == ML_ERR_SINGULAR,
                     "zero polynomial divisor is rejected");

        ml_polyq_clear (&x);
        ml_polyq_clear (&f);
        ml_polyq_clear (&g);
        ml_polyq_clear (&quotient);
        ml_polyq_clear (&remainder);
        ml_polyq_clear (&product);
        ml_polyq_clear (&derivative);
        ml_polyq_clear (&gcd);
        ml_mpq_clear (&at_two);
        ml_mpq_clear (&value);
    }

    return ml_test_summary (&ctx);
}
