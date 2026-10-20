#include <ml_math_extra.h>
#include <ml_exact.h>
#include <ml_polyq.h>
#include <ml_gfp.h>
#include <ml_school.h>
int main (void) {
    int quotient = 0;
    double remainder = ml_remquo (0x1p100, 3.0, &quotient);
    ml_mpz_t n, d;
    ml_mpq_t q, x, value, expected;
    ml_interval_t enclosure;
    ml_polyq_t p;
    ml_gfp_ctx_t field;
    ml_school_solution_set_t school_set;
    uint32_t inverse = 0u;
    double metres = 0.0;
    int failed;

    ml_mpz_init (&n);
    ml_mpz_init (&d);
    ml_mpq_init (&q);
    ml_mpq_init (&x);
    ml_mpq_init (&value);
    ml_mpq_init (&expected);
    ml_polyq_init (&p);
    ml_gfp_ctx_init (&field);
    school_set.count = 0u;
    failed = !(remainder == 1.0 && quotient == 85);
    failed |= ml_mpq_set_decimal (&expected, "0.125") != ML_SUCCESS;
    failed |= ml_mpq_to_interval (&expected, &enclosure) != ML_SUCCESS ||
              enclosure.lo != 0.125 || enclosure.hi != 0.125;
    failed |= ml_mpz_set_u64 (&n, 1u) != ML_SUCCESS;
    failed |= ml_mpz_set_u64 (&d, 3u) != ML_SUCCESS;
    failed |= ml_mpq_set_ratio (&q, &n, &d) != ML_SUCCESS;
    failed |= ml_mpq_to_interval (&q, &enclosure) != ML_SUCCESS;
    failed |= !ml_interval_contains (enclosure, 1.0 / 3.0);
    failed |= ml_mpq_set_i64 (&x, 2) != ML_SUCCESS;
    failed |= ml_polyq_set_coefficient (&p, 1u, &q) != ML_SUCCESS;
    failed |= ml_polyq_evaluate (&value, &p, &x) != ML_SUCCESS;
    failed |= ml_mpz_set_u64 (&n, 2u) != ML_SUCCESS;
    failed |= ml_mpq_set_ratio (&expected, &n, &d) != ML_SUCCESS;
    failed |= ml_mpq_compare (&value, &expected, &quotient) != ML_SUCCESS;
    failed |= quotient != 0;
    failed |= ml_gfp_ctx_set_prime (&field, 7u) != ML_SUCCESS;
    failed |= ml_gfp_inv (&field, 3u, &inverse) != ML_SUCCESS || inverse != 5u;
    failed |= ml_school_convert_length (1.0, ML_SCHOOL_LENGTH_FT,
                                        ML_SCHOOL_LENGTH_M, &metres) != ML_SUCCESS ||
              metres != 0.3048;
    failed |= ml_school_solve_quadratic_inequality (1.0, 0.0, -1.0,
                                                    ML_SCHOOL_REL_LE, &school_set) != ML_SUCCESS ||
              school_set.count != 1u ||
              school_set.interval [0].lower != -1.0 || school_set.interval [0].upper != 1.0;

    ml_polyq_clear (&p);
    ml_mpq_clear (&q);
    ml_mpq_clear (&x);
    ml_mpq_clear (&value);
    ml_mpq_clear (&expected);
    ml_mpz_clear (&n);
    ml_mpz_clear (&d);
    return failed;
}
