#include "ml_school.h"

#include "ml_core.h"
#include <float.h>
#include <limits.h>
#include <math.h>

#define ML_SCHOOL_MAX_EXACT_INDEX UINT64_C (9007199254740992)
static ml_status_t school_store (long double value, double *out) {
    double converted;
    if (!out) return ML_ERR_INVALID_ARG;
    if (!(value <= (long double) DBL_MAX &&
          value >= -(long double) DBL_MAX)) return ML_ERR_OVERFLOW;
    converted = (double) value;
    if (!ml_isfinite (converted)) return ML_ERR_OVERFLOW;
    *out = converted;
    return ML_SUCCESS;
}
static ml_status_t school_convert (double value, int from, int to,
                                   int count, const double *to_base,
                                   double *out) {
    long double converted;
    if (!out || !to_base || from < 0 || to < 0 || from >= count || to >= count)
        return ML_ERR_INVALID_ARG;
    if (!ml_isfinite (value)) return ML_ERR_NAN_INPUT;
    if (value < 0.0) return ML_ERR_INVALID_ARG;
    converted = (long double) value * (long double) to_base [from] /
                (long double) to_base [to];
    return school_store (converted, out);
}
ml_status_t ml_school_convert_length (double value,
                                      ml_school_length_unit_t from, ml_school_length_unit_t to, double *out) {
    static const double metres [] = {
        0.001, 0.01, 1.0, 1000.0, 0.0254, 0.3048, 0.9144, 1609.344};
    return school_convert (value, (int) from, (int) to,
                           (int) (sizeof (metres) / sizeof (metres [0])), metres, out);
}
ml_status_t ml_school_convert_mass (double value,
                                    ml_school_mass_unit_t from, ml_school_mass_unit_t to, double *out) {
    static const double kilograms [] = {
        1e-6, 0.001, 1.0, 1000.0, 0.028349523125, 0.45359237};
    return school_convert (value, (int) from, (int) to,
                           (int) (sizeof (kilograms) / sizeof (kilograms [0])), kilograms, out);
}
ml_status_t ml_school_convert_time (double value,
                                    ml_school_time_unit_t from, ml_school_time_unit_t to, double *out) {
    static const double seconds [] = {
        0.001, 1.0, 60.0, 3600.0, 86400.0, 604800.0, 31536000.0};
    return school_convert (value, (int) from, (int) to,
                           (int) (sizeof (seconds) / sizeof (seconds [0])), seconds, out);
}
ml_status_t ml_school_convert_volume (double value,
                                      ml_school_volume_unit_t from, ml_school_volume_unit_t to, double *out) {
    /* Values are millilitres; customary liquid measures are US units. */
    static const double millilitres [] = {
        1.0, 1000.0, 29.5735295625, 236.5882365,
        473.176473, 946.352946, 3785.411784,
        0.001, 1.0, 1000000.0, 16.387064, 28316.846592};
    return school_convert (value, (int) from, (int) to,
                           (int) (sizeof (millilitres) / sizeof (millilitres [0])), millilitres, out);
}
ml_status_t ml_school_convert_area (double value,
                                    ml_school_area_unit_t from, ml_school_area_unit_t to, double *out) {
    static const double square_metres [] = {
        1e-6, 1e-4, 1.0, 1e6, 10000.0,
        0.00064516, 0.09290304, 4046.8564224, 2589988.110336};
    return school_convert (value, (int) from, (int) to,
                           (int) (sizeof (square_metres) / sizeof (square_metres [0])),
                           square_metres, out);
}
ml_status_t ml_school_convert_temperature (double value,
                                           ml_school_temperature_unit_t from, ml_school_temperature_unit_t to,
                                           double *out) {
    long double celsius;
    long double converted;
    if (!out || (int) from < 0 || (int) from > (int) ML_SCHOOL_TEMP_K ||
        (int) to < 0 || (int) to > (int) ML_SCHOOL_TEMP_K)
        return ML_ERR_INVALID_ARG;
    if (!ml_isfinite (value)) return ML_ERR_NAN_INPUT;
    if ((from == ML_SCHOOL_TEMP_C && value < -273.15) ||
        (from == ML_SCHOOL_TEMP_F && value < -459.67) ||
        (from == ML_SCHOOL_TEMP_K && value < 0.0)) return ML_ERR_INVALID_ARG;
    if (from == ML_SCHOOL_TEMP_C) celsius = (long double) value;
    else if (from == ML_SCHOOL_TEMP_F) {
        if (value == -459.67) celsius = -273.15L;
        else celsius = ((long double) value - 32.0L) * (5.0L / 9.0L);
    } else {
        if (value == 0.0) celsius = -273.15L;
        else celsius = (long double) value - 273.15L;
    }
    if (to == ML_SCHOOL_TEMP_C) converted = celsius;
    else if (to == ML_SCHOOL_TEMP_F)
        converted = celsius * (9.0L / 5.0L) + 32.0L;
    else converted = celsius + 273.15L;
    return school_store (converted, out);
}
ml_status_t ml_school_percent_of (double value, double percent, double *out) {
    if (!out) return ML_ERR_INVALID_ARG;
    if (!ml_isfinite (value) || !ml_isfinite (percent)) return ML_ERR_NAN_INPUT;
    return school_store ((long double) value * (long double) percent / 100.0L, out);
}
ml_status_t ml_school_percent_change (double old_value, double new_value,
                                      double *percent_change) {
    long double result;
    if (!percent_change) return ML_ERR_INVALID_ARG;
    if (!ml_isfinite (old_value) || !ml_isfinite (new_value))
        return ML_ERR_NAN_INPUT;
    if (old_value == 0.0) return ML_ERR_SINGULAR;
    result = 100.0L * ((long double) new_value - (long double) old_value) /
             (long double) old_value;
    return school_store (result, percent_change);
}
static int school_round_mode_valid (ml_school_round_mode_t mode) {
    return (int) mode >= (int) ML_SCHOOL_ROUND_HALF_AWAY_FROM_ZERO &&
           (int) mode <= (int) ML_SCHOOL_ROUND_CEILING;
}
static long double school_round_value (long double value,
                                       ml_school_round_mode_t mode) {
    long double lower, fraction;
    switch (mode) {
    case ML_SCHOOL_ROUND_HALF_AWAY_FROM_ZERO:
        return copysignl (floorl (fabsl (value) + 0.5L), value);
    case ML_SCHOOL_ROUND_HALF_TO_EVEN:
        lower = floorl (value);
        fraction = value - lower;
        if (fraction < 0.5L) return lower;
        if (fraction > 0.5L) return lower + 1.0L;
        return fmodl (fabsl (lower), 2.0L) == 0.0L ? lower : lower + 1.0L;
    case ML_SCHOOL_ROUND_TOWARD_ZERO:
        return truncl (value);
    case ML_SCHOOL_ROUND_FLOOR:
        return floorl (value);
    case ML_SCHOOL_ROUND_CEILING:
        return ceill (value);
    default:
        return value;
    }
}
ml_status_t ml_school_round_decimal (double value, int decimal_places,
                                     ml_school_round_mode_t mode, double *out) {
    long double scale, scaled, rounded;
    if (!out || !school_round_mode_valid (mode)) return ML_ERR_INVALID_ARG;
    if (!ml_isfinite (value)) return ML_ERR_NAN_INPUT;
    if (decimal_places < -308 || decimal_places > 308)
        return ML_ERR_INVALID_ARG;
    if (value == 0.0) {
        *out = value;
        return ML_SUCCESS;
    }
    scale = powl (10.0L, (long double) decimal_places);
    if (scale == 0.0L) {
        *out = 0.0;
        return ML_SUCCESS;
    }
    if (!isfinite (scale)) {
        /* Finer than the input's representable spacing: rounding is a no-op. */
        *out = value;
        return ML_SUCCESS;
    }
    scaled = (long double) value * scale;
    if (!isfinite (scaled)) {
        if (decimal_places > 0) {
            *out = value;
            return ML_SUCCESS;
        }
        return ML_ERR_OVERFLOW;
    }
    rounded = school_round_value (scaled, mode);
    return school_store (rounded / scale, out);
}
ml_status_t ml_school_round_significant (double value, int digits,
                                         ml_school_round_mode_t mode, double *out) {
    double exponent_value;
    int exponent, decimal_places;
    if (!out || digits < 1 || digits > 17 || !school_round_mode_valid (mode))
        return ML_ERR_INVALID_ARG;
    if (!ml_isfinite (value)) return ML_ERR_NAN_INPUT;
    if (value == 0.0) {
        *out = value;
        return ML_SUCCESS;
    }
    exponent_value = floor (log10 (fabs (value)));
    if (!ml_isfinite (exponent_value) || exponent_value < (double) INT_MIN ||
        exponent_value > (double) INT_MAX) return ML_ERR_UNSUPPORTED;
    exponent = (int) exponent_value;
    decimal_places = digits - 1 - exponent;
    if (decimal_places < -308 || decimal_places > 308)
        return ML_ERR_UNSUPPORTED;
    return ml_school_round_decimal (value, decimal_places, mode, out);
}
static int school_index_ok (uint64_t index) {
    return index <= ML_SCHOOL_MAX_EXACT_INDEX;
}
ml_status_t ml_school_arithmetic_term (double first, double difference,
                                       uint64_t index, double *out) {
    if (!out) return ML_ERR_INVALID_ARG;
    if (!ml_isfinite (first) || !ml_isfinite (difference)) return ML_ERR_NAN_INPUT;
    if (!school_index_ok (index)) return ML_ERR_UNSUPPORTED;
    return school_store ((long double) first +
                             (long double) difference * (long double) index,
                         out);
}
ml_status_t ml_school_arithmetic_sum (double first, double difference,
                                      uint64_t terms, double *out) {
    long double count;
    long double last;
    if (!out) return ML_ERR_INVALID_ARG;
    if (!ml_isfinite (first) || !ml_isfinite (difference)) return ML_ERR_NAN_INPUT;
    if (!school_index_ok (terms)) return ML_ERR_UNSUPPORTED;
    if (terms == 0u) return school_store (0.0L, out);
    count = (long double) terms;
    last = (long double) first + (long double) difference * (count - 1.0L);
    return school_store (count * ((long double) first + last) / 2.0L, out);
}
ml_status_t ml_school_geometric_term (double first, double ratio,
                                      uint64_t index, double *out) {
    double power;
    if (!out) return ML_ERR_INVALID_ARG;
    if (!ml_isfinite (first) || !ml_isfinite (ratio)) return ML_ERR_NAN_INPUT;
    if (!school_index_ok (index)) return ML_ERR_UNSUPPORTED;
    if (first == 0.0) {
        *out = first;
        return ML_SUCCESS;
    }
    power = pow (ratio, (double) index);
    if (!ml_isfinite (power)) return ML_ERR_OVERFLOW;
    return school_store ((long double) first * (long double) power, out);
}
ml_status_t ml_school_geometric_sum (double first, double ratio,
                                     uint64_t terms, double *out) {
    uint64_t remaining = terms;
    double accumulated_power = 1.0;
    double accumulated_sum = 0.0;
    double block_power = ratio;
    double block_sum = 1.0;
    if (!out) return ML_ERR_INVALID_ARG;
    if (!ml_isfinite (first) || !ml_isfinite (ratio)) return ML_ERR_NAN_INPUT;
    if (!school_index_ok (terms)) return ML_ERR_UNSUPPORTED;
    if (terms == 0u || first == 0.0) {
        *out = 0.0;
        return ML_SUCCESS;
    }
    while (remaining != 0u) {
        if ((remaining & UINT64_C (1)) != 0u) {
            accumulated_sum += accumulated_power * block_sum;
            accumulated_power *= block_power;
            if (!ml_isfinite (accumulated_sum) || !ml_isfinite (accumulated_power))
                return ML_ERR_OVERFLOW;
        }
        remaining >>= 1;
        if (remaining != 0u) {
            block_sum *= 1.0 + block_power;
            block_power *= block_power;
            if (!ml_isfinite (block_sum) || !ml_isfinite (block_power))
                return ML_ERR_OVERFLOW;
        }
    }
    return school_store ((long double) first * (long double) accumulated_sum, out);
}
ml_status_t ml_school_simple_amount (double principal, double rate_per_period,
                                     uint64_t periods, double *amount) {
    long double result;
    if (!amount) return ML_ERR_INVALID_ARG;
    if (!ml_isfinite (principal) || !ml_isfinite (rate_per_period))
        return ML_ERR_NAN_INPUT;
    if (principal < 0.0) return ML_ERR_INVALID_ARG;
    if (!school_index_ok (periods)) return ML_ERR_UNSUPPORTED;
    result = (long double) principal *
             (1.0L + (long double) rate_per_period * (long double) periods);
    return school_store (result, amount);
}
ml_status_t ml_school_compound_amount (double principal, double rate_per_period,
                                       uint64_t periods, double *amount) {
    double factor;
    if (!amount) return ML_ERR_INVALID_ARG;
    if (!ml_isfinite (principal) || !ml_isfinite (rate_per_period))
        return ML_ERR_NAN_INPUT;
    if (principal < 0.0 || rate_per_period < -1.0) return ML_ERR_INVALID_ARG;
    if (!school_index_ok (periods)) return ML_ERR_UNSUPPORTED;
    if (principal == 0.0) {
        *amount = 0.0;
        return ML_SUCCESS;
    }
    factor = pow (1.0 + rate_per_period, (double) periods);
    if (!ml_isfinite (factor)) return ML_ERR_OVERFLOW;
    return school_store ((long double) principal * (long double) factor, amount);
}
ml_status_t ml_school_solve_linear (double a, double b, double *x) {
    if (!x) return ML_ERR_INVALID_ARG;
    if (!ml_isfinite (a) || !ml_isfinite (b)) return ML_ERR_NAN_INPUT;
    if (a == 0.0) return ML_ERR_SINGULAR;
    return school_store (-(long double) b / (long double) a, x);
}
ml_status_t ml_school_slope (double x1, double y1, double x2, double y2,
                             double *slope) {
    long double result;
    if (!slope) return ML_ERR_INVALID_ARG;
    if (!ml_isfinite (x1) || !ml_isfinite (y1) ||
        !ml_isfinite (x2) || !ml_isfinite (y2)) return ML_ERR_NAN_INPUT;
    if (x1 == x2) return ML_ERR_SINGULAR;
    result = ((long double) y2 - (long double) y1) /
             ((long double) x2 - (long double) x1);
    return school_store (result, slope);
}
ml_status_t ml_school_line_intersection (double a1, double b1, double c1,
                                         double a2, double b2, double c2, double *x, double *y) {
    long double s1, s2, la1, lb1, lc1, la2, lb2, lc2;
    long double determinant, result_x, result_y;
    double out_x, out_y;
    if (!x || !y || x == y) return ML_ERR_INVALID_ARG;
    if (!ml_isfinite (a1) || !ml_isfinite (b1) || !ml_isfinite (c1) ||
        !ml_isfinite (a2) || !ml_isfinite (b2) || !ml_isfinite (c2))
        return ML_ERR_NAN_INPUT;
    s1 = fmax (fabs (a1), fabs (b1));
    s2 = fmax (fabs (a2), fabs (b2));
    if (s1 == 0.0L || s2 == 0.0L) return ML_ERR_INVALID_ARG;
    la1 = (long double) a1 / s1;
    lb1 = (long double) b1 / s1;
    lc1 = (long double) c1 / s1;
    la2 = (long double) a2 / s2;
    lb2 = (long double) b2 / s2;
    lc2 = (long double) c2 / s2;
    determinant = la1 * lb2 - la2 * lb1;
    if (determinant == 0.0L) return ML_ERR_SINGULAR;
    result_x = (lc1 * lb2 - lc2 * lb1) / determinant;
    result_y = (la1 * lc2 - la2 * lc1) / determinant;
    if (school_store (result_x, &out_x) != ML_SUCCESS ||
        school_store (result_y, &out_y) != ML_SUCCESS) return ML_ERR_OVERFLOW;
    *x = out_x;
    *y = out_y;
    return ML_SUCCESS;
}
static int school_relation_valid (ml_school_relation_t relation) {
    return (int) relation >= (int) ML_SCHOOL_REL_LT &&
           (int) relation <= (int) ML_SCHOOL_REL_GE;
}
static int school_relation_accepts (int sign, ml_school_relation_t relation) {
    switch (relation) {
    case ML_SCHOOL_REL_LT:
        return sign < 0;
    case ML_SCHOOL_REL_LE:
        return sign <= 0;
    case ML_SCHOOL_REL_GT:
        return sign > 0;
    case ML_SCHOOL_REL_GE:
        return sign >= 0;
    default:
        return 0;
    }
}
static void school_set_empty (ml_school_solution_set_t *set) {
    set->count = 0u;
    set->interval [0].lower = 0.0;
    set->interval [0].upper = 0.0;
    set->interval [0].lower_closed = 0;
    set->interval [0].upper_closed = 0;
    set->interval [1] = set->interval [0];
}
static void school_set_all (ml_school_solution_set_t *set) {
    school_set_empty (set);
    set->count = 1u;
    set->interval [0].lower = -HUGE_VAL;
    set->interval [0].upper = HUGE_VAL;
}
static ml_status_t school_store_root (long double root, double *out) {
    return school_store (root, out);
}
ml_status_t ml_school_solve_linear_inequality (double a, double b,
                                               ml_school_relation_t relation, ml_school_solution_set_t *out) {
    ml_school_solution_set_t result;
    double root;
    int sign, accepts_zero;
    if (!out || !school_relation_valid (relation)) return ML_ERR_INVALID_ARG;
    if (!ml_isfinite (a) || !ml_isfinite (b)) return ML_ERR_NAN_INPUT;
    school_set_empty (&result);
    if (a == 0.0) {
        sign = (b > 0.0) ? 1 : ((b < 0.0) ? -1 : 0);
        if (school_relation_accepts (sign, relation)) school_set_all (&result);
        *out = result;
        return ML_SUCCESS;
    }
    if (school_store_root (-(long double) b / (long double) a, &root) != ML_SUCCESS)
        return ML_ERR_OVERFLOW;
    accepts_zero = school_relation_accepts (0, relation);
    result.count = 1u;
    if ((a > 0.0 && relation <= ML_SCHOOL_REL_LE) ||
        (a < 0.0 && relation >= ML_SCHOOL_REL_GE)) {
        result.interval [0].lower = -HUGE_VAL;
        result.interval [0].upper = root;
        result.interval [0].lower_closed = 0;
        result.interval [0].upper_closed = accepts_zero;
    } else {
        result.interval [0].lower = root;
        result.interval [0].upper = HUGE_VAL;
        result.interval [0].lower_closed = accepts_zero;
        result.interval [0].upper_closed = 0;
    }
    *out = result;
    return ML_SUCCESS;
}
static ml_status_t school_quadratic_root (long double a, long double b,
                                          long double c, long double discriminant, int upper_root, double *out) {
    long double sqrt_discriminant = sqrtl (discriminant);
    long double q = -0.5L * (b + copysignl (sqrt_discriminant, b));
    long double root;
    if (q == 0.0L) {
        root = -b / (2.0L * a);
    } else if (upper_root) {
        root = c / q;
    } else {
        root = q / a;
    }
    return school_store_root (root, out);
}
ml_status_t ml_school_solve_quadratic_inequality (double a, double b,
                                                  double c, ml_school_relation_t relation, ml_school_solution_set_t *out) {
    ml_school_solution_set_t result;
    long double scale, la, lb, lc, discriminant;
    int outer_sign, accepts_outer, accepts_inner, accepts_zero;
    double root1, root2;
    if (!out || !school_relation_valid (relation)) return ML_ERR_INVALID_ARG;
    if (!ml_isfinite (a) || !ml_isfinite (b) || !ml_isfinite (c))
        return ML_ERR_NAN_INPUT;
    if (a == 0.0) return ml_school_solve_linear_inequality (b, c, relation, out);
    scale = (long double) fmax (fabs (a), fmax (fabs (b), fabs (c)));
    la = (long double) a / scale;
    lb = (long double) b / scale;
    lc = (long double) c / scale;
    if (la == 0.0L || (b != 0.0 && lb == 0.0L) ||
        (c != 0.0 && lc == 0.0L)) return ML_ERR_UNSUPPORTED;
    discriminant = fmal (lb, lb, -4.0L * la * lc);
    school_set_empty (&result);
    outer_sign = (a > 0.0) ? 1 : -1;
    accepts_outer = school_relation_accepts (outer_sign, relation);
    accepts_inner = school_relation_accepts (-outer_sign, relation);
    accepts_zero = school_relation_accepts (0, relation);
    if (discriminant < 0.0L) {
        if (accepts_outer) school_set_all (&result);
        *out = result;
        return ML_SUCCESS;
    }
    if (!isfinite (discriminant)) return ML_ERR_OVERFLOW;
    if (discriminant == 0.0L) {
        if (school_quadratic_root (la, lb, lc, 0.0L, 0, &root1) != ML_SUCCESS)
            return ML_ERR_OVERFLOW;
        if (accepts_outer && accepts_zero) {
            school_set_all (&result);
        } else if (accepts_outer) {
            result.count = 2u;
            result.interval [0].lower = -HUGE_VAL;
            result.interval [0].upper = root1;
            result.interval [0].lower_closed = 0;
            result.interval [0].upper_closed = 0;
            result.interval [1].lower = root1;
            result.interval [1].upper = HUGE_VAL;
            result.interval [1].lower_closed = 0;
            result.interval [1].upper_closed = 0;
        } else if (accepts_zero) {
            result.count = 1u;
            result.interval [0].lower = root1;
            result.interval [0].upper = root1;
            result.interval [0].lower_closed = 1;
            result.interval [0].upper_closed = 1;
        }
        *out = result;
        return ML_SUCCESS;
    }
    if (school_quadratic_root (la, lb, lc, discriminant, 0, &root1) != ML_SUCCESS ||
        school_quadratic_root (la, lb, lc, discriminant, 1, &root2) != ML_SUCCESS)
        return ML_ERR_OVERFLOW;
    if (root1 > root2) {
        double swap = root1;
        root1 = root2;
        root2 = swap;
    }
    if (accepts_outer) {
        result.count = 2u;
        result.interval [0].lower = -HUGE_VAL;
        result.interval [0].upper = root1;
        result.interval [0].lower_closed = 0;
        result.interval [0].upper_closed = accepts_zero;
        result.interval [1].lower = root2;
        result.interval [1].upper = HUGE_VAL;
        result.interval [1].lower_closed = accepts_zero;
        result.interval [1].upper_closed = 0;
    } else if (accepts_inner) {
        result.count = 1u;
        result.interval [0].lower = root1;
        result.interval [0].upper = root2;
        result.interval [0].lower_closed = accepts_zero;
        result.interval [0].upper_closed = accepts_zero;
    }
    *out = result;
    return ML_SUCCESS;
}
int ml_school_solution_set_contains (const ml_school_solution_set_t *set,
                                     double x) {
    size_t i;
    if (!set || !ml_isfinite (x) || set->count > 2u) return 0;
    for (i = 0u; i < set->count; ++i) {
        const ml_school_interval_t *interval = &set->interval [i];
        int lower_ok, upper_ok;
        if (ml_isnan (interval->lower) || ml_isnan (interval->upper) ||
            interval->lower > interval->upper ||
            (interval->lower == -HUGE_VAL && interval->lower_closed) ||
            (interval->upper == HUGE_VAL && interval->upper_closed)) return 0;
        lower_ok = (x > interval->lower) ||
                   (interval->lower_closed && x == interval->lower);
        upper_ok = (x < interval->upper) ||
                   (interval->upper_closed && x == interval->upper);
        if (lower_ok && upper_ok) return 1;
    }
    return 0;
}
static ml_status_t school_nonnegative3 (double a, double b, double c) {
    if (!ml_isfinite (a) || !ml_isfinite (b) || !ml_isfinite (c))
        return ML_ERR_NAN_INPUT;
    if (a < 0.0 || b < 0.0 || c < 0.0) return ML_ERR_INVALID_ARG;
    return ML_SUCCESS;
}
ml_status_t ml_school_rectangle_area (double width, double height, double *area) {
    ml_status_t st;
    if (!area) return ML_ERR_INVALID_ARG;
    st = school_nonnegative3 (width, height, 0.0);
    if (st != ML_SUCCESS) return st;
    return school_store ((long double) width * (long double) height, area);
}
ml_status_t ml_school_rectangle_perimeter (double width, double height,
                                           double *perimeter) {
    ml_status_t st;
    if (!perimeter) return ML_ERR_INVALID_ARG;
    st = school_nonnegative3 (width, height, 0.0);
    if (st != ML_SUCCESS) return st;
    return school_store (2.0L * ((long double) width + (long double) height),
                         perimeter);
}
ml_status_t ml_school_triangle_area (double base, double height, double *area) {
    ml_status_t st;
    if (!area) return ML_ERR_INVALID_ARG;
    st = school_nonnegative3 (base, height, 0.0);
    if (st != ML_SUCCESS) return st;
    return school_store ((long double) base * (long double) height / 2.0L, area);
}
ml_status_t ml_school_triangle_perimeter (double side_a, double side_b,
                                          double side_c, double *perimeter) {
    double largest;
    long double other_sum;
    ml_status_t st;
    if (!perimeter) return ML_ERR_INVALID_ARG;
    st = school_nonnegative3 (side_a, side_b, side_c);
    if (st != ML_SUCCESS) return st;
    largest = fmax (side_a, fmax (side_b, side_c));
    if (largest == 0.0) return ML_ERR_INVALID_ARG;
    if (largest == side_a)
        other_sum = (long double) side_b + (long double) side_c;
    else if (largest == side_b)
        other_sum = (long double) side_a + (long double) side_c;
    else
        other_sum = (long double) side_a + (long double) side_b;
    if (!(other_sum > (long double) largest)) return ML_ERR_INVALID_ARG;
    return school_store ((long double) side_a + (long double) side_b +
                             (long double) side_c,
                         perimeter);
}
ml_status_t ml_school_trapezoid_area (double parallel_a, double parallel_b,
                                      double height, double *area) {
    ml_status_t st;
    if (!area) return ML_ERR_INVALID_ARG;
    st = school_nonnegative3 (parallel_a, parallel_b, height);
    if (st != ML_SUCCESS) return st;
    return school_store (((long double) parallel_a + (long double) parallel_b) *
                             (long double) height / 2.0L,
                         area);
}
ml_status_t ml_school_circle_area (double radius, double *area) {
    ml_status_t st;
    if (!area) return ML_ERR_INVALID_ARG;
    st = school_nonnegative3 (radius, 0.0, 0.0);
    if (st != ML_SUCCESS) return st;
    return school_store ((long double) ML_PI * (long double) radius *
                             (long double) radius,
                         area);
}
ml_status_t ml_school_circle_circumference (double radius,
                                            double *circumference) {
    ml_status_t st;
    if (!circumference) return ML_ERR_INVALID_ARG;
    st = school_nonnegative3 (radius, 0.0, 0.0);
    if (st != ML_SUCCESS) return st;
    return school_store (2.0L * (long double) ML_PI * (long double) radius,
                         circumference);
}
ml_status_t ml_school_circle_sector_area (double radius, double angle_radians,
                                          double *area) {
    ml_status_t st;
    if (!area) return ML_ERR_INVALID_ARG;
    st = school_nonnegative3 (radius, angle_radians, 0.0);
    if (st != ML_SUCCESS) return st;
    if ((long double) angle_radians > 2.0L * (long double) ML_PI)
        return ML_ERR_INVALID_ARG;
    return school_store (0.5L * (long double) radius * (long double) radius *
                             (long double) angle_radians,
                         area);
}
ml_status_t ml_school_rectangular_prism_volume (double length, double width,
                                                double height, double *volume) {
    ml_status_t st;
    if (!volume) return ML_ERR_INVALID_ARG;
    st = school_nonnegative3 (length, width, height);
    if (st != ML_SUCCESS) return st;
    return school_store ((long double) length * (long double) width *
                             (long double) height,
                         volume);
}
ml_status_t ml_school_rectangular_prism_surface_area (double length,
                                                      double width, double height, double *area) {
    ml_status_t st;
    long double l, w, h;
    if (!area) return ML_ERR_INVALID_ARG;
    st = school_nonnegative3 (length, width, height);
    if (st != ML_SUCCESS) return st;
    l = (long double) length;
    w = (long double) width;
    h = (long double) height;
    return school_store (2.0L * (l * w + w * h + h * l), area);
}
ml_status_t ml_school_pyramid_volume (double base_area, double height,
                                      double *volume) {
    ml_status_t st;
    if (!volume) return ML_ERR_INVALID_ARG;
    st = school_nonnegative3 (base_area, height, 0.0);
    if (st != ML_SUCCESS) return st;
    return school_store ((long double) base_area * (long double) height / 3.0L,
                         volume);
}
ml_status_t ml_school_cylinder_surface_area (double radius, double height,
                                             double *area) {
    ml_status_t st;
    if (!area) return ML_ERR_INVALID_ARG;
    st = school_nonnegative3 (radius, height, 0.0);
    if (st != ML_SUCCESS) return st;
    return school_store (2.0L * (long double) ML_PI * (long double) radius *
                             ((long double) radius + (long double) height),
                         area);
}
