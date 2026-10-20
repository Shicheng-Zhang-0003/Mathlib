#include "test_harness.h"
#include "ml_exact.h"
#include "ml_interval.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
static int mpz_equals_string (const ml_mpz_t *a, const char *expected) {
    char buf [1600];
    int n = ml_mpz_to_string (a, buf, sizeof (buf), 10u);
    return n >= 0 && strcmp (buf, expected) == 0;
}
int main (void) {
    ml_test_ctx_t ctx;
    ml_test_init (&ctx, "Exact and certified foundations");

    {
        ml_mpz_t a, b, sum, product, q, rem, ten, rebuilt, g;
        ml_mpz_init (&a);
        ml_mpz_init (&b);
        ml_mpz_init (&sum);
        ml_mpz_init (&product);
        ml_mpz_init (&q);
        ml_mpz_init (&rem);
        ml_mpz_init (&ten);
        ml_mpz_init (&rebuilt);
        ml_mpz_init (&g);
        ASSERT_TRUE (&ctx, ml_mpz_set_string (&a, "123456789012345678901234567890", 10u) == ML_SUCCESS,
                     "parse a value beyond native integer range");
        ASSERT_TRUE (&ctx, ml_mpz_set_string (&b, "9876543210", 10u) == ML_SUCCESS,
                     "parse second exact integer");
        ASSERT_TRUE (&ctx, ml_mpz_add (&sum, &a, &b) == ML_SUCCESS && mpz_equals_string (&sum, "123456789012345678911111111100"),
                     "arbitrary-precision addition");
        ASSERT_TRUE (&ctx, ml_mpz_mul (&product, &a, &b) == ML_SUCCESS,
                     "arbitrary-precision multiplication");
        ASSERT_TRUE (&ctx, ml_mpz_divrem (&q, &rem, &product, &a) == ML_SUCCESS && ml_mpz_compare (&q, &b) == 0 && ml_mpz_is_zero (&rem),
                     "product division identity");
        ASSERT_TRUE (&ctx, ml_mpz_set_i64 (&a, -100) == ML_SUCCESS && ml_mpz_set_u64 (&b, 7u) == ML_SUCCESS && ml_mpz_divrem (&q, &rem, &a, &b) == ML_SUCCESS && mpz_equals_string (&q, "-14") && mpz_equals_string (&rem, "-2"),
                     "signed division truncates toward zero");
        ASSERT_TRUE (&ctx, ml_mpz_set_u64 (&a, 48u) == ML_SUCCESS && ml_mpz_set_u64 (&b, 18u) == ML_SUCCESS && ml_mpz_gcd (&g, &a, &b) == ML_SUCCESS && mpz_equals_string (&g, "6"), "arbitrary-precision gcd");
        {
            uint64_t state = UINT64_C (0x9e3779b97f4a7c15);
            int i;
            for (i = 0; i < 100; i++) {
                int64_t x, y, expected;
                char expected_text [64];
                state ^= state >> 12;
                state ^= state << 25;
                state ^= state >> 27;
                x = (int64_t) ((state * UINT64_C (0x2545f4914f6cdd1d)) % UINT64_C (2000000001)) -
                    INT64_C (1000000000);
                state ^= state >> 12;
                state ^= state << 25;
                state ^= state >> 27;
                y = (int64_t) ((state * UINT64_C (0x2545f4914f6cdd1d)) % UINT64_C (2000000001)) -
                    INT64_C (1000000000);
                (void) ml_mpz_set_i64 (&a, x);
                (void) ml_mpz_set_i64 (&b, y);
                expected = x + y;
                (void) snprintf (expected_text, sizeof (expected_text), "%lld",
                                 (long long) expected);
                ASSERT_TRUE (&ctx, ml_mpz_add (&sum, &a, &b) == ML_SUCCESS && mpz_equals_string (&sum, expected_text),
                             "random exact addition agrees with int64 oracle");
                expected = x * y;
                (void) snprintf (expected_text, sizeof (expected_text), "%lld",
                                 (long long) expected);
                ASSERT_TRUE (&ctx, ml_mpz_mul (&product, &a, &b) == ML_SUCCESS && mpz_equals_string (&product, expected_text),
                             "random exact multiplication agrees with int64 oracle");
                if (y != 0) {
                    expected = x / y;
                    (void) snprintf (expected_text, sizeof (expected_text), "%lld",
                                     (long long) expected);
                    ASSERT_TRUE (&ctx, ml_mpz_divrem (&q, &rem, &a, &b) == ML_SUCCESS && mpz_equals_string (&q, expected_text),
                                 "random exact quotient agrees with int64 oracle");
                    expected = x % y;
                    (void) snprintf (expected_text, sizeof (expected_text), "%lld",
                                     (long long) expected);
                    ASSERT_TRUE (&ctx, mpz_equals_string (&rem, expected_text),
                                 "random exact remainder agrees with int64 oracle");
                }
            }
        }
        ASSERT_TRUE (&ctx, ml_mpz_set_u64 (&ten, 10u) == ML_SUCCESS,
                     "prepare decimal divisor");
        ml_mpz_clear (&a);
        ml_mpz_init (&a);
        {
            char *digits = (char *) malloc (1501u);
            char *expected = (char *) malloc (1500u);
            size_t i;
            if (digits == NULL || expected == NULL) {
                ASSERT_TRUE (&ctx, 0, "allocate long-integer test input");
            } else {
                memset (digits, '9', 1500u);
                digits [1500] = '\0';
                memset (expected, '9', 1499u);
                expected [1499] = '\0';
                ASSERT_TRUE (&ctx, ml_mpz_set_string (&a, digits, 10u) == ML_SUCCESS && ml_mpz_bit_length (&a) > 4096u,
                             "integer grows past legacy 4096-bit capacity");
                ASSERT_TRUE (&ctx, ml_mpz_divrem (&q, &rem, &a, &ten) == ML_SUCCESS && mpz_equals_string (&rem, "9"),
                             "large-integer division remainder");
                {
                    char *buf = (char *) malloc (1500u);
                    int n = buf ? ml_mpz_to_string (&q, buf, 1500u, 10u) : -1;
                    ASSERT_TRUE (&ctx, n == 1499 && buf != NULL && strcmp (buf, expected) == 0,
                                 "large quotient has every expected digit");
                    free (buf);
                }
                ASSERT_TRUE (&ctx, ml_mpz_mul (&rebuilt, &q, &ten) == ML_SUCCESS && ml_mpz_add (&rebuilt, &rebuilt, &rem) == ML_SUCCESS && ml_mpz_compare (&rebuilt, &a) == 0,
                             "large division reconstructs dividend");
                for (i = 0u; i < 1500u; i++) {
                    if (digits [i] != '9') break;
                }
                ASSERT_TRUE (&ctx, i == 1500u,
                             "test string remains intact after parsing");
                ASSERT_TRUE (&ctx, ml_mpz_mul (&product, &a, &a) == ML_SUCCESS && ml_mpz_divrem (&q, &rem, &product, &a) == ML_SUCCESS && ml_mpz_compare (&q, &a) == 0 && ml_mpz_is_zero (&rem),
                             "large multi-limb product and quotient identity");
            }
            free (digits);
            free (expected);
        }
        {
            char before [1600], after [1600];
            (void) ml_mpz_set_u64 (&a, 314u);
            (void) ml_mpz_to_string (&a, before, sizeof (before), 10u);
            ASSERT_TRUE (&ctx, ml_mpz_set_string (&a, "12z", 10u) == ML_ERR_INVALID_ARG && ml_mpz_to_string (&a, after, sizeof (after), 10u) == 3 && strcmp (before, after) == 0,
                         "invalid parse leaves destination unchanged");
        }
        ml_mpz_clear (&a);
        ml_mpz_clear (&b);
        ml_mpz_clear (&sum);
        ml_mpz_clear (&product);
        ml_mpz_clear (&q);
        ml_mpz_clear (&rem);
        ml_mpz_clear (&ten);
        ml_mpz_clear (&rebuilt);
        ml_mpz_clear (&g);
    }

    {
        ml_mpz_t n, d;
        ml_mpq_t a, b, c, expected;
        int cmp = 99;
        ml_mpz_init (&n);
        ml_mpz_init (&d);
        ml_mpq_init (&a);
        ml_mpq_init (&b);
        ml_mpq_init (&c);
        ml_mpq_init (&expected);
        (void) ml_mpz_set_u64 (&n, 100u);
        (void) ml_mpz_set_u64 (&d, 250u);
        ASSERT_TRUE (&ctx, ml_mpq_set_ratio (&a, &n, &d) == ML_SUCCESS && mpz_equals_string (&a.numerator, "2") && mpz_equals_string (&a.denominator, "5"),
                     "rational normalization reduces common factors");
        ASSERT_TRUE (&ctx, ml_mpq_set_decimal (&a, "0.125") == ML_SUCCESS && mpz_equals_string (&a.numerator, "1") && mpz_equals_string (&a.denominator, "8"),
                     "exact decimal parses to a reduced fraction");
        ASSERT_TRUE (&ctx, ml_mpq_set_decimal (&a, "-12.50e-1") == ML_SUCCESS && mpz_equals_string (&a.numerator, "-5") && mpz_equals_string (&a.denominator, "4"),
                     "signed decimal with exponent parses exactly");
        ASSERT_TRUE (&ctx, ml_mpq_set_decimal (&a, ".375") == ML_SUCCESS && mpz_equals_string (&a.numerator, "3") && mpz_equals_string (&a.denominator, "8"),
                     "decimal may omit the leading zero");
        ASSERT_TRUE (&ctx, ml_mpq_set_decimal (&a, "1.") == ML_SUCCESS && mpz_equals_string (&a.numerator, "1") && mpz_equals_string (&a.denominator, "1"),
                     "decimal may omit fractional digits");
        ASSERT_TRUE (&ctx, ml_mpq_set_decimal (&a, "1e3") == ML_SUCCESS && mpz_equals_string (&a.numerator, "1000") && mpz_equals_string (&a.denominator, "1"),
                     "positive exponent shifts decimal scale exactly");
        (void) ml_mpq_set_i64 (&a, 7);
        ASSERT_TRUE (&ctx, ml_mpq_set_decimal (&a, "1.2.3") == ML_ERR_INVALID_ARG && mpz_equals_string (&a.numerator, "7") && mpz_equals_string (&a.denominator, "1"),
                     "malformed decimal leaves rational output unchanged");
        ASSERT_TRUE (&ctx, ml_mpq_set_decimal (&a, "1e10001") == ML_ERR_UNSUPPORTED && mpz_equals_string (&a.numerator, "7"),
                     "decimal exponent resource cap is explicit");
        (void) ml_mpz_set_u64 (&n, 1u);
        (void) ml_mpz_set_u64 (&d, 3u);
        ASSERT_TRUE (&ctx, ml_mpq_set_ratio (&a, &n, &d) == ML_SUCCESS,
                     "construct one third");
        (void) ml_mpz_set_u64 (&n, 1u);
        (void) ml_mpz_set_u64 (&d, 6u);
        ASSERT_TRUE (&ctx, ml_mpq_set_ratio (&b, &n, &d) == ML_SUCCESS,
                     "construct one sixth");
        ASSERT_TRUE (&ctx, ml_mpq_add (&c, &a, &b) == ML_SUCCESS && mpz_equals_string (&c.numerator, "1") && mpz_equals_string (&c.denominator, "2"),
                     "exact rational addition canonicalizes");
        ASSERT_TRUE (&ctx, ml_mpq_mul (&c, &a, &b) == ML_SUCCESS && mpz_equals_string (&c.numerator, "1") && mpz_equals_string (&c.denominator, "18"),
                     "exact rational multiplication");
        ASSERT_TRUE (&ctx, ml_mpq_div (&c, &a, &b) == ML_SUCCESS && mpz_equals_string (&c.numerator, "2") && mpz_equals_string (&c.denominator, "1"),
                     "exact rational division");
        (void) ml_mpz_set_u64 (&n, 1u);
        (void) ml_mpz_set_u64 (&d, 2u);
        (void) ml_mpq_set_ratio (&expected, &n, &d);
        ASSERT_TRUE (&ctx, ml_mpq_compare (&a, &expected, &cmp) == ML_SUCCESS && cmp < 0,
                     "exact rational comparison");
        ASSERT_TRUE (&ctx, ml_mpq_add (&a, &a, &b) == ML_SUCCESS && mpz_equals_string (&a.numerator, "1") && mpz_equals_string (&a.denominator, "2"),
                     "rational arithmetic supports output aliasing");
        (void) ml_mpz_set_u64 (&d, 0u);
        ASSERT_TRUE (&ctx, ml_mpq_set_ratio (&c, &n, &d) == ML_ERR_SINGULAR && mpz_equals_string (&c.numerator, "2") && mpz_equals_string (&c.denominator, "1"),
                     "zero denominator rejected without changing output");
        if (ml_interval_is_supported ()) {
            ml_interval_t enclosure = {-1.0, 1.0};
            ASSERT_TRUE (&ctx, ml_mpq_to_interval (&b, &enclosure) == ML_SUCCESS && ml_interval_contains (enclosure, 1.0 / 6.0) && enclosure.lo < enclosure.hi && nextafter (enclosure.lo, INFINITY) == enclosure.hi,
                         "rational converts to tight adjacent-double enclosure");
        } else {
            ASSERT_TRUE (&ctx, 1, "rational enclosure support is explicitly optional");
        }
        ml_mpz_clear (&n);
        ml_mpz_clear (&d);
        ml_mpq_clear (&a);
        ml_mpq_clear (&b);
        ml_mpq_clear (&c);
        ml_mpq_clear (&expected);
    }

    {
        ml_interval_t a, b, r, keep;
        ASSERT_TRUE (&ctx, ml_interval_is_supported (),
                     "runtime reports IEEE binary64 round-to-nearest support");
        ASSERT_TRUE (&ctx, ml_interval_make (0.1, 0.2, &a) == ML_SUCCESS && ml_interval_point (0.2, &b) == ML_SUCCESS && ml_interval_add (a, b, &r) == ML_SUCCESS && ml_interval_contains (r, 0.3),
                     "outward interval addition contains exact decimal sum");
        ASSERT_TRUE (&ctx, ml_interval_make (-1.0, 2.0, &a) == ML_SUCCESS && ml_interval_make (3.0, 4.0, &b) == ML_SUCCESS && ml_interval_mul (a, b, &r) == ML_SUCCESS && r.lo <= -3.0 && r.hi >= 8.0,
                     "interval multiplication encloses all endpoint products");
        ASSERT_TRUE (&ctx, ml_interval_point (1.0, &a) == ML_SUCCESS && ml_interval_point (3.0, &b) == ML_SUCCESS && ml_interval_div (a, b, &r) == ML_SUCCESS && ml_interval_contains (r, 1.0 / 3.0),
                     "interval division encloses one third");
        (void) ml_interval_make (-1.0, 1.0, &b);
        keep.lo = -7.0;
        keep.hi = 9.0;
        r = keep;
        ASSERT_TRUE (&ctx, ml_interval_div (a, b, &r) == ML_ERR_SINGULAR && r.lo == keep.lo && r.hi == keep.hi,
                     "division by interval containing zero fails atomically");
        (void) ml_interval_point (DBL_MAX, &a);
        (void) ml_interval_point (2.0, &b);
        r = keep;
        ASSERT_TRUE (&ctx, ml_interval_mul (a, b, &r) == ML_ERR_OVERFLOW && r.lo == keep.lo && r.hi == keep.hi,
                     "unrepresentable finite enclosure fails atomically");
        ASSERT_TRUE (&ctx, ml_interval_make (2.0, 1.0, &r) == ML_ERR_INVALID_ARG && r.lo == keep.lo && r.hi == keep.hi,
                     "reversed interval rejected without modifying output");
    }

    return ml_test_summary (&ctx);
}
