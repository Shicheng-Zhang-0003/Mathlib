#ifndef MATHLIB_ML_SCHOOL_H
#define MATHLIB_ML_SCHOOL_H

#include "ml_compiler.h"
#include "ml_types.h"
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif
/* Practical elementary-through-high-school mathematics.
 *
 * Conversions and formula routines accept finite inputs and write their
 * output only on success. Geometric dimensions, masses, times and volumes
 * must be nonnegative. Temperature conversion treats its input as an
 * absolute temperature and rejects values below absolute zero.
 *
 * School-model routines use binary64 arithmetic; they are numerical
 * conveniences, not exact decimal/financial accounting. Use ml_mpq_t for
 * exact rational calculations. */

typedef enum {
    ML_SCHOOL_LENGTH_MM = 0,
    ML_SCHOOL_LENGTH_CM,
    ML_SCHOOL_LENGTH_M,
    ML_SCHOOL_LENGTH_KM,
    ML_SCHOOL_LENGTH_IN,
    ML_SCHOOL_LENGTH_FT,
    ML_SCHOOL_LENGTH_YD,
    ML_SCHOOL_LENGTH_MI
} ml_school_length_unit_t;
typedef enum {
    ML_SCHOOL_MASS_MG = 0,
    ML_SCHOOL_MASS_G,
    ML_SCHOOL_MASS_KG,
    ML_SCHOOL_MASS_TONNE,
    ML_SCHOOL_MASS_OZ,
    ML_SCHOOL_MASS_LB
} ml_school_mass_unit_t;
typedef enum {
    ML_SCHOOL_TIME_MS = 0,
    ML_SCHOOL_TIME_S,
    ML_SCHOOL_TIME_MIN,
    ML_SCHOOL_TIME_HOUR,
    ML_SCHOOL_TIME_DAY,
    ML_SCHOOL_TIME_WEEK,
    ML_SCHOOL_TIME_YEAR_365D
} ml_school_time_unit_t;
typedef enum {
    ML_SCHOOL_VOLUME_ML = 0,
    ML_SCHOOL_VOLUME_L,
    ML_SCHOOL_VOLUME_US_FL_OZ,
    ML_SCHOOL_VOLUME_US_CUP,
    ML_SCHOOL_VOLUME_US_PINT,
    ML_SCHOOL_VOLUME_US_QUART,
    ML_SCHOOL_VOLUME_US_GALLON,
    ML_SCHOOL_VOLUME_MM3,
    ML_SCHOOL_VOLUME_CM3,
    ML_SCHOOL_VOLUME_M3,
    ML_SCHOOL_VOLUME_IN3,
    ML_SCHOOL_VOLUME_FT3
} ml_school_volume_unit_t;
typedef enum {
    ML_SCHOOL_AREA_MM2 = 0,
    ML_SCHOOL_AREA_CM2,
    ML_SCHOOL_AREA_M2,
    ML_SCHOOL_AREA_KM2,
    ML_SCHOOL_AREA_HECTARE,
    ML_SCHOOL_AREA_IN2,
    ML_SCHOOL_AREA_FT2,
    ML_SCHOOL_AREA_ACRE,
    ML_SCHOOL_AREA_MI2
} ml_school_area_unit_t;
typedef enum {
    ML_SCHOOL_TEMP_C = 0,
    ML_SCHOOL_TEMP_F,
    ML_SCHOOL_TEMP_K
} ml_school_temperature_unit_t;
ML_API ml_status_t ml_school_convert_length (double value,
                                             ml_school_length_unit_t from, ml_school_length_unit_t to, double *out);
ML_API ml_status_t ml_school_convert_mass (double value,
                                           ml_school_mass_unit_t from, ml_school_mass_unit_t to, double *out);
ML_API ml_status_t ml_school_convert_time (double value,
                                           ml_school_time_unit_t from, ml_school_time_unit_t to, double *out);
ML_API ml_status_t ml_school_convert_volume (double value,
                                             ml_school_volume_unit_t from, ml_school_volume_unit_t to, double *out);
ML_API ml_status_t ml_school_convert_area (double value,
                                           ml_school_area_unit_t from, ml_school_area_unit_t to, double *out);
ML_API ml_status_t ml_school_convert_temperature (double value,
                                                  ml_school_temperature_unit_t from, ml_school_temperature_unit_t to,
                                                  double *out);

/* Percent change is 100*(new-old)/old; old == 0 has no finite percentage. */
ML_API ml_status_t ml_school_percent_of (double value, double percent,
                                         double *out);
ML_API ml_status_t ml_school_percent_change (double old_value,
                                             double new_value, double *percent_change);
typedef enum {
    ML_SCHOOL_ROUND_HALF_AWAY_FROM_ZERO = 0,
    ML_SCHOOL_ROUND_HALF_TO_EVEN,
    ML_SCHOOL_ROUND_TOWARD_ZERO,
    ML_SCHOOL_ROUND_FLOOR,
    ML_SCHOOL_ROUND_CEILING
} ml_school_round_mode_t;
/* Rounds the exact binary64 input after scaling by 10^decimal_places, with
 * an explicit tie rule. decimal_places is limited to [-308,308]. */
ML_API ml_status_t ml_school_round_decimal (double value, int decimal_places,
                                            ml_school_round_mode_t mode, double *out);
/* Significant digits are limited to [1,17]. Returns UNSUPPORTED if the
 * equivalent decimal-place count falls outside [-308,308]. Rounding is
 * applied to the represented binary64 value, so familiar decimal literals
 * such as 1.005 retain their binary representation behavior. */
ML_API ml_status_t ml_school_round_significant (double value, int digits,
                                                ml_school_round_mode_t mode, double *out);

/* Sequences are indexed from zero. `terms` counts a0 through a_(terms-1);
 * zero terms sum to 0. Indices/counts above 2^53 are rejected because a
 * binary64 result cannot reliably distinguish all such integer indices. */
ML_API ml_status_t ml_school_arithmetic_term (double first, double difference,
                                              uint64_t index, double *out);
ML_API ml_status_t ml_school_arithmetic_sum (double first, double difference,
                                             uint64_t terms, double *out);
ML_API ml_status_t ml_school_geometric_term (double first, double ratio,
                                             uint64_t index, double *out);
ML_API ml_status_t ml_school_geometric_sum (double first, double ratio,
                                            uint64_t terms, double *out);

/* Simple/compound growth use a per-period rate: 0.05 means 5% each period.
 * Compound rate must be >= -1.0. These formula helpers are not an accounting
 * or tax library and do not model fees, inflation, or calendar conventions. */
ML_API ml_status_t ml_school_simple_amount (double principal,
                                            double rate_per_period, uint64_t periods, double *amount);
ML_API ml_status_t ml_school_compound_amount (double principal,
                                              double rate_per_period, uint64_t periods, double *amount);

/* Linear equation a*x + b = 0. A non-unique or inconsistent equation has no
 * single answer and returns ML_ERR_SINGULAR. */
ML_API ml_status_t ml_school_solve_linear (double a, double b, double *x);
ML_API ml_status_t ml_school_slope (double x1, double y1, double x2, double y2,
                                    double *slope);
/* Intersect lines A*x + B*y = C. Parallel/coincident lines return
 * ML_ERR_SINGULAR. */
ML_API ml_status_t ml_school_line_intersection (double a1, double b1,
                                                double c1, double a2, double b2, double c2, double *x, double *y);
typedef enum {
    ML_SCHOOL_REL_LT = 0,
    ML_SCHOOL_REL_LE,
    ML_SCHOOL_REL_GT,
    ML_SCHOOL_REL_GE
} ml_school_relation_t;
typedef struct {
    double lower;
    double upper;
    int lower_closed;
    int upper_closed;
} ml_school_interval_t;
typedef struct {
    size_t count; /* 0: empty, 1: one interval, 2: two intervals */
    ml_school_interval_t interval [2];
} ml_school_solution_set_t;
/* Solve a*x+b relation 0 and a*x^2+b*x+c relation 0. Infinite interval
 * endpoints are represented by +/-INFINITY and are always open. Quadratic
 * roots are binary64 approximations, not certified enclosures; a discriminant
 * close to zero can be classified differently by different FP platforms. */
ML_API ml_status_t ml_school_solve_linear_inequality (double a, double b,
                                                      ml_school_relation_t relation, ml_school_solution_set_t *out);
ML_API ml_status_t ml_school_solve_quadratic_inequality (double a, double b,
                                                         double c, ml_school_relation_t relation, ml_school_solution_set_t *out);
ML_API int ml_school_solution_set_contains (const ml_school_solution_set_t *set,
                                            double x);

ML_API ml_status_t ml_school_rectangle_area (double width, double height,
                                             double *area);
ML_API ml_status_t ml_school_rectangle_perimeter (double width, double height,
                                                  double *perimeter);
ML_API ml_status_t ml_school_triangle_area (double base, double height,
                                            double *area);
ML_API ml_status_t ml_school_triangle_perimeter (double side_a, double side_b,
                                                 double side_c, double *perimeter);
ML_API ml_status_t ml_school_trapezoid_area (double parallel_a,
                                             double parallel_b, double height, double *area);
ML_API ml_status_t ml_school_circle_area (double radius, double *area);
ML_API ml_status_t ml_school_circle_circumference (double radius,
                                                   double *circumference);
/* Sector angle is in radians and must lie in [0, 2*pi]. */
ML_API ml_status_t ml_school_circle_sector_area (double radius,
                                                 double angle_radians, double *area);
ML_API ml_status_t ml_school_rectangular_prism_volume (double length,
                                                       double width, double height, double *volume);
ML_API ml_status_t ml_school_rectangular_prism_surface_area (double length,
                                                             double width, double height, double *area);
ML_API ml_status_t ml_school_pyramid_volume (double base_area, double height,
                                             double *volume);
ML_API ml_status_t ml_school_cylinder_surface_area (double radius,
                                                    double height, double *area);

#ifdef __cplusplus
}
#endif
#endif /* MATHLIB_ML_SCHOOL_H */
