#include "test_harness.h"
#include "ml_school.h"

#include <math.h>
#include <stdint.h>
int main (void) {
    ml_test_ctx_t ctx;
    double value = 123.0, x = 123.0, y = 456.0;
    ml_test_init (&ctx, "Pre-university mathematics utilities");

    ASSERT_TRUE (&ctx, ml_school_convert_length (5.0, ML_SCHOOL_LENGTH_FT, ML_SCHOOL_LENGTH_M, &value) == ML_SUCCESS,
                 "convert feet to metres");
    ASSERT_NEAR (&ctx, value, 1.524, 1e-15, "five feet in metres");
    ASSERT_TRUE (&ctx, ml_school_convert_volume (1.0, ML_SCHOOL_VOLUME_US_GALLON, ML_SCHOOL_VOLUME_L, &value) == ML_SUCCESS,
                 "convert US gallon to litre");
    ASSERT_NEAR (&ctx, value, 3.785411784, 1e-14, "US liquid gallon in litres");
    ASSERT_TRUE (&ctx, ml_school_convert_volume (1000.0, ML_SCHOOL_VOLUME_CM3, ML_SCHOOL_VOLUME_L, &value) == ML_SUCCESS,
                 "convert cubic centimetres to litres");
    ASSERT_NEAR (&ctx, value, 1.0, 0.0, "1000 cubic centimetres in one litre");
    ASSERT_TRUE (&ctx, ml_school_convert_volume (1.0, ML_SCHOOL_VOLUME_FT3, ML_SCHOOL_VOLUME_L, &value) == ML_SUCCESS,
                 "convert cubic feet to litres");
    ASSERT_NEAR (&ctx, value, 28.316846592, 1e-13,
                 "cubic foot in litres");
    ASSERT_TRUE (&ctx, ml_school_convert_area (1.0, ML_SCHOOL_AREA_HECTARE, ML_SCHOOL_AREA_M2, &value) == ML_SUCCESS,
                 "convert hectares to square metres");
    ASSERT_NEAR (&ctx, value, 10000.0, 0.0, "hectare area conversion");
    ASSERT_TRUE (&ctx, ml_school_convert_area (1.0, ML_SCHOOL_AREA_ACRE, ML_SCHOOL_AREA_FT2, &value) == ML_SUCCESS,
                 "convert acres to square feet");
    ASSERT_NEAR (&ctx, value, 43560.0, 1e-10, "acre area conversion");
    ASSERT_TRUE (&ctx, ml_school_convert_mass (16.0, ML_SCHOOL_MASS_OZ, ML_SCHOOL_MASS_LB, &value) == ML_SUCCESS,
                 "convert ounces to pounds");
    ASSERT_NEAR (&ctx, value, 1.0, 0.0, "sixteen ounces in one pound");
    ASSERT_TRUE (&ctx, ml_school_convert_time (2.0, ML_SCHOOL_TIME_WEEK, ML_SCHOOL_TIME_DAY, &value) == ML_SUCCESS,
                 "convert weeks to days");
    ASSERT_NEAR (&ctx, value, 14.0, 0.0, "two weeks in days");
    ASSERT_TRUE (&ctx, ml_school_convert_temperature (32.0, ML_SCHOOL_TEMP_F, ML_SCHOOL_TEMP_C, &value) == ML_SUCCESS,
                 "convert freezing point to Celsius");
    ASSERT_NEAR (&ctx, value, 0.0, 1e-14, "Fahrenheit freezing point");
    ASSERT_TRUE (&ctx, ml_school_convert_temperature (-459.67, ML_SCHOOL_TEMP_F, ML_SCHOOL_TEMP_K, &value) == ML_SUCCESS,
                 "convert Fahrenheit absolute zero");
    ASSERT_NEAR (&ctx, value, 0.0, 0.0, "Fahrenheit absolute zero is 0 K");
    value = 77.0;
    ASSERT_TRUE (&ctx, ml_school_convert_length (-1.0, ML_SCHOOL_LENGTH_M, ML_SCHOOL_LENGTH_CM, &value) == ML_ERR_INVALID_ARG && value == 77.0,
                 "negative length rejected atomically");
    ASSERT_TRUE (&ctx, ml_school_convert_temperature (-300.0, ML_SCHOOL_TEMP_C, ML_SCHOOL_TEMP_K, &value) == ML_ERR_INVALID_ARG && value == 77.0,
                 "temperature below absolute zero rejected atomically");
    ASSERT_TRUE (&ctx, ml_school_convert_time (1.0, (ml_school_time_unit_t) 99, ML_SCHOOL_TIME_S, &value) == ML_ERR_INVALID_ARG && value == 77.0,
                 "invalid unit rejected atomically");

    ASSERT_TRUE (&ctx, ml_school_percent_of (80.0, 12.5, &value) == ML_SUCCESS,
                 "calculate a percentage of a quantity");
    ASSERT_NEAR (&ctx, value, 10.0, 0.0, "12.5 percent of 80");
    ASSERT_TRUE (&ctx, ml_school_percent_change (80.0, 100.0, &value) == ML_SUCCESS,
                 "calculate percentage change");
    ASSERT_NEAR (&ctx, value, 25.0, 0.0, "increase from 80 to 100");
    value = 9.0;
    ASSERT_TRUE (&ctx, ml_school_percent_change (0.0, 4.0, &value) == ML_ERR_SINGULAR && value == 9.0,
                 "percentage change from zero has no finite result");

    ASSERT_TRUE (&ctx, ml_school_round_decimal (2.5, 0, ML_SCHOOL_ROUND_HALF_AWAY_FROM_ZERO, &value) == ML_SUCCESS,
                 "round half away from zero");
    ASSERT_NEAR (&ctx, value, 3.0, 0.0, "positive half rounds away from zero");
    ASSERT_TRUE (&ctx, ml_school_round_decimal (-2.5, 0, ML_SCHOOL_ROUND_HALF_AWAY_FROM_ZERO, &value) == ML_SUCCESS,
                 "round negative half away from zero");
    ASSERT_NEAR (&ctx, value, -3.0, 0.0, "negative half rounds away from zero");
    ASSERT_TRUE (&ctx, ml_school_round_decimal (2.5, 0, ML_SCHOOL_ROUND_HALF_TO_EVEN, &value) == ML_SUCCESS,
                 "round ties to even");
    ASSERT_NEAR (&ctx, value, 2.0, 0.0, "2.5 rounds to the even integer 2");
    ASSERT_TRUE (&ctx, ml_school_round_decimal (3.5, 0, ML_SCHOOL_ROUND_HALF_TO_EVEN, &value) == ML_SUCCESS,
                 "round odd tie to even");
    ASSERT_NEAR (&ctx, value, 4.0, 0.0, "3.5 rounds to the even integer 4");
    ASSERT_TRUE (&ctx, ml_school_round_decimal (1.25, 1, ML_SCHOOL_ROUND_HALF_TO_EVEN, &value) == ML_SUCCESS,
                 "round decimal place with explicit tie mode");
    ASSERT_NEAR (&ctx, value, 1.2, 1e-15, "1.25 to one decimal ties to even");
    ASSERT_TRUE (&ctx, ml_school_round_decimal (-2.9, 0, ML_SCHOOL_ROUND_TOWARD_ZERO, &value) == ML_SUCCESS,
                 "round toward zero");
    ASSERT_NEAR (&ctx, value, -2.0, 0.0, "toward-zero rounding");
    ASSERT_TRUE (&ctx, ml_school_round_significant (12345.0, 3, ML_SCHOOL_ROUND_HALF_AWAY_FROM_ZERO, &value) == ML_SUCCESS,
                 "round to significant digits");
    ASSERT_NEAR (&ctx, value, 12300.0, 0.0, "three significant figures");
    ASSERT_TRUE (&ctx, ml_school_round_significant (9.99, 2, ML_SCHOOL_ROUND_HALF_AWAY_FROM_ZERO, &value) == ML_SUCCESS,
                 "significant rounding carries into next decade");
    ASSERT_NEAR (&ctx, value, 10.0, 0.0, "rounding 9.99 to two significant digits");
    value = 9.0;
    ASSERT_TRUE (&ctx, ml_school_round_decimal (1.0, 0, (ml_school_round_mode_t) 99, &value) == ML_ERR_INVALID_ARG && value == 9.0, "invalid rounding mode rejected atomically");

    ASSERT_TRUE (&ctx, ml_school_arithmetic_term (5.0, 3.0, 4u, &value) == ML_SUCCESS, "arithmetic sequence term");
    ASSERT_NEAR (&ctx, value, 17.0, 0.0, "a4 for 5, 8, 11, ...");
    ASSERT_TRUE (&ctx, ml_school_arithmetic_sum (5.0, 3.0, 5u, &value) == ML_SUCCESS, "arithmetic sequence sum");
    ASSERT_NEAR (&ctx, value, 55.0, 0.0, "sum of first five arithmetic terms");
    ASSERT_TRUE (&ctx, ml_school_arithmetic_sum (1.0, 1.0, 0u, &value) == ML_SUCCESS && value == 0.0, "empty sequence sum is zero");
    ASSERT_TRUE (&ctx, ml_school_geometric_term (2.0, 3.0, 3u, &value) == ML_SUCCESS, "geometric sequence term");
    ASSERT_NEAR (&ctx, value, 54.0, 0.0, "a3 for 2, 6, 18, 54");
    ASSERT_TRUE (&ctx, ml_school_geometric_sum (2.0, 3.0, 4u, &value) == ML_SUCCESS, "geometric sequence sum");
    ASSERT_NEAR (&ctx, value, 80.0, 0.0, "sum of first four geometric terms");
    ASSERT_TRUE (&ctx, ml_school_geometric_sum (7.0, 1.0, 5u, &value) == ML_SUCCESS, "geometric sum handles ratio one");
    ASSERT_NEAR (&ctx, value, 35.0, 0.0, "constant geometric sequence sum");
    ASSERT_TRUE (&ctx, ml_school_geometric_sum (2.0, -1.0, 6u, &value) == ML_SUCCESS && value == 0.0, "alternating geometric sum");
    value = 9.0;
    ASSERT_TRUE (&ctx, ml_school_geometric_term (1.0, 2.0, UINT64_C (9007199254740993), &value) == ML_ERR_UNSUPPORTED && value == 9.0, "inexact binary64 index rejected atomically");

    ASSERT_TRUE (&ctx, ml_school_simple_amount (1000.0, 0.05, 4u, &value) == ML_SUCCESS, "simple interest amount");
    ASSERT_NEAR (&ctx, value, 1200.0, 1e-12, "simple interest formula");
    ASSERT_TRUE (&ctx, ml_school_compound_amount (1000.0, 0.05, 2u, &value) == ML_SUCCESS, "compound interest amount");
    ASSERT_NEAR (&ctx, value, 1102.5, 1e-12, "compound interest formula");
    value = 9.0;
    ASSERT_TRUE (&ctx, ml_school_compound_amount (10.0, -1.1, 2u, &value) == ML_ERR_INVALID_ARG && value == 9.0,
                 "compound rate below total loss rejected atomically");

    ASSERT_TRUE (&ctx, ml_school_solve_linear (2.0, 6.0, &x) == ML_SUCCESS,
                 "solve one-variable linear equation");
    ASSERT_NEAR (&ctx, x, -3.0, 0.0, "solution of 2x+6=0");
    {
        ml_school_solution_set_t solutions;
        ASSERT_TRUE (&ctx, ml_school_solve_linear_inequality (2.0, 4.0, ML_SCHOOL_REL_LT, &solutions) == ML_SUCCESS && solutions.count == 1u && solutions.interval [0].lower == -HUGE_VAL && solutions.interval [0].upper == -2.0 && !solutions.interval [0].upper_closed,
                     "strict linear inequality has an open endpoint");
        ASSERT_TRUE (&ctx, ml_school_solution_set_contains (&solutions, -3.0) && !ml_school_solution_set_contains (&solutions, -2.0) && !ml_school_solution_set_contains (&solutions, 0.0),
                     "linear inequality interval membership");
        ASSERT_TRUE (&ctx, ml_school_solve_linear_inequality (-2.0, 4.0, ML_SCHOOL_REL_GE, &solutions) == ML_SUCCESS && solutions.interval [0].upper == 2.0 && solutions.interval [0].upper_closed,
                     "inequality direction reverses for negative coefficient");
        ASSERT_TRUE (&ctx, ml_school_solve_linear_inequality (0.0, 1.0, ML_SCHOOL_REL_LT, &solutions) == ML_SUCCESS && solutions.count == 0u,
                     "false constant inequality has empty solution set");
        ASSERT_TRUE (&ctx, ml_school_solve_linear_inequality (0.0, 0.0, ML_SCHOOL_REL_LE, &solutions) == ML_SUCCESS && solutions.count == 1u && ml_school_solution_set_contains (&solutions, -1e100) && ml_school_solution_set_contains (&solutions, 1e100),
                     "true constant inequality covers all real numbers");

        ASSERT_TRUE (&ctx, ml_school_solve_quadratic_inequality (1.0, 0.0, -1.0, ML_SCHOOL_REL_LT, &solutions) == ML_SUCCESS && solutions.count == 1u && solutions.interval [0].lower == -1.0 && solutions.interval [0].upper == 1.0 && !solutions.interval [0].lower_closed && !solutions.interval [0].upper_closed,
                     "quadratic strict-between-roots solution");
        ASSERT_TRUE (&ctx, !ml_school_solution_set_contains (&solutions, -1.0) && ml_school_solution_set_contains (&solutions, 0.0) && !ml_school_solution_set_contains (&solutions, 1.0),
                     "strict quadratic roots are excluded");
        ASSERT_TRUE (&ctx, ml_school_solve_quadratic_inequality (1.0, 0.0, -1.0, ML_SCHOOL_REL_GE, &solutions) == ML_SUCCESS && solutions.count == 2u && solutions.interval [0].upper_closed && solutions.interval [1].lower_closed,
                     "quadratic outside-roots solution includes roots");
        ASSERT_TRUE (&ctx, ml_school_solution_set_contains (&solutions, -2.0) && ml_school_solution_set_contains (&solutions, -1.0) && !ml_school_solution_set_contains (&solutions, 0.0) && ml_school_solution_set_contains (&solutions, 1.0) && ml_school_solution_set_contains (&solutions, 2.0),
                     "two-interval quadratic membership");
        ASSERT_TRUE (&ctx, ml_school_solve_quadratic_inequality (1.0, 0.0, 1.0, ML_SCHOOL_REL_GT, &solutions) == ML_SUCCESS && solutions.count == 1u && ml_school_solution_set_contains (&solutions, 0.0),
                     "positive quadratic with no real roots is always positive");
        ASSERT_TRUE (&ctx, ml_school_solve_quadratic_inequality (1.0, 0.0, 1.0, ML_SCHOOL_REL_LE, &solutions) == ML_SUCCESS && solutions.count == 0u,
                     "positive quadratic cannot be nonpositive");
        ASSERT_TRUE (&ctx, ml_school_solve_quadratic_inequality (1.0, -2.0, 1.0, ML_SCHOOL_REL_GT, &solutions) == ML_SUCCESS && solutions.count == 2u && !ml_school_solution_set_contains (&solutions, 1.0),
                     "strict inequality excludes a repeated root");
        ASSERT_TRUE (&ctx, ml_school_solve_quadratic_inequality (1.0, -2.0, 1.0, ML_SCHOOL_REL_LE, &solutions) == ML_SUCCESS && solutions.count == 1u && solutions.interval [0].lower == 1.0 && solutions.interval [0].upper == 1.0 && ml_school_solution_set_contains (&solutions, 1.0),
                     "non-strict inequality reduces to repeated-root singleton");
        ASSERT_TRUE (&ctx, ml_school_solve_quadratic_inequality (1.0, -1e16, 1.0, ML_SCHOOL_REL_LE, &solutions) == ML_SUCCESS && solutions.count == 1u && solutions.interval [0].lower > 0.0 && solutions.interval [0].lower < 1e-15 && solutions.interval [0].upper >= 1e16 && ml_school_solution_set_contains (&solutions, 1.0),
                     "quadratic root calculation avoids cancellation");
        ASSERT_TRUE (&ctx, ml_school_solve_quadratic_inequality (0.0, 2.0, -4.0, ML_SCHOOL_REL_LT, &solutions) == ML_SUCCESS && solutions.interval [0].upper == 2.0,
                     "quadratic solver handles linear degeneration");
        solutions.count = 77u;
        ASSERT_TRUE (&ctx, ml_school_solve_quadratic_inequality (1.0, 2.0, 1.0, (ml_school_relation_t) 99, &solutions) == ML_ERR_INVALID_ARG && solutions.count == 77u,
                     "invalid relation rejected without replacing output");
    }
    ASSERT_TRUE (&ctx, ml_school_slope (1.0, 2.0, 5.0, 10.0, &value) == ML_SUCCESS, "slope between two points");
    ASSERT_NEAR (&ctx, value, 2.0, 0.0, "slope of line through points");
    value = 9.0;
    ASSERT_TRUE (&ctx, ml_school_slope (2.0, 1.0, 2.0, 4.0, &value) == ML_ERR_SINGULAR && value == 9.0,
                 "vertical line has undefined slope");
    ASSERT_TRUE (&ctx, ml_school_line_intersection (1.0, 1.0, 3.0, 1.0, -1.0, 1.0, &x, &y) == ML_SUCCESS,
                 "intersect two lines");
    ASSERT_NEAR (&ctx, x, 2.0, 0.0, "line intersection x coordinate");
    ASSERT_NEAR (&ctx, y, 1.0, 0.0, "line intersection y coordinate");
    x = 9.0;
    y = 8.0;
    ASSERT_TRUE (&ctx, ml_school_line_intersection (1.0, 2.0, 3.0, 2.0, 4.0, 6.0, &x, &y) == ML_ERR_SINGULAR && x == 9.0 && y == 8.0, "coincident lines have no unique intersection");

    ASSERT_TRUE (&ctx, ml_school_rectangle_area (3.0, 4.0, &value) == ML_SUCCESS,
                 "rectangle area");
    ASSERT_NEAR (&ctx, value, 12.0, 0.0, "rectangle area formula");
    ASSERT_TRUE (&ctx, ml_school_rectangle_perimeter (3.0, 4.0, &value) == ML_SUCCESS, "rectangle perimeter");
    ASSERT_NEAR (&ctx, value, 14.0, 0.0, "rectangle perimeter formula");
    ASSERT_TRUE (&ctx, ml_school_triangle_area (6.0, 4.0, &value) == ML_SUCCESS,
                 "triangle area from base and height");
    ASSERT_NEAR (&ctx, value, 12.0, 0.0, "triangle area formula");
    ASSERT_TRUE (&ctx, ml_school_triangle_perimeter (3.0, 4.0, 5.0, &value) == ML_SUCCESS, "valid triangle perimeter");
    ASSERT_NEAR (&ctx, value, 12.0, 0.0, "triangle perimeter formula");
    value = 9.0;
    ASSERT_TRUE (&ctx, ml_school_triangle_perimeter (1.0, 2.0, 3.0, &value) == ML_ERR_INVALID_ARG && value == 9.0,
                 "degenerate triangle rejected atomically");
    ASSERT_TRUE (&ctx, ml_school_trapezoid_area (4.0, 8.0, 3.0, &value) == ML_SUCCESS, "trapezoid area");
    ASSERT_NEAR (&ctx, value, 18.0, 0.0, "trapezoid area formula");
    ASSERT_TRUE (&ctx, ml_school_circle_area (2.0, &value) == ML_SUCCESS,
                 "circle area");
    ASSERT_NEAR (&ctx, value, 4.0 * ML_PI, 0.0, "circle area formula");
    ASSERT_TRUE (&ctx, ml_school_circle_circumference (2.0, &value) == ML_SUCCESS,
                 "circle circumference");
    ASSERT_NEAR (&ctx, value, 4.0 * ML_PI, 0.0, "circumference formula");
    ASSERT_TRUE (&ctx, ml_school_circle_sector_area (2.0, ML_PI / 2.0, &value) == ML_SUCCESS, "circle sector area");
    ASSERT_NEAR (&ctx, value, ML_PI, 0.0, "quarter-circle sector area");
    ASSERT_TRUE (&ctx, ml_school_rectangular_prism_volume (2.0, 3.0, 4.0, &value) == ML_SUCCESS, "rectangular prism volume");
    ASSERT_NEAR (&ctx, value, 24.0, 0.0, "prism volume formula");
    ASSERT_TRUE (&ctx, ml_school_rectangular_prism_surface_area (2.0, 3.0, 4.0, &value) == ML_SUCCESS, "rectangular prism surface area");
    ASSERT_NEAR (&ctx, value, 52.0, 0.0, "prism surface area formula");
    ASSERT_TRUE (&ctx, ml_school_pyramid_volume (12.0, 5.0, &value) == ML_SUCCESS, "pyramid volume");
    ASSERT_NEAR (&ctx, value, 20.0, 0.0, "pyramid volume formula");
    ASSERT_TRUE (&ctx, ml_school_cylinder_surface_area (2.0, 3.0, &value) == ML_SUCCESS, "cylinder total surface area");
    ASSERT_NEAR (&ctx, value, 20.0 * ML_PI, 0.0, "cylinder surface formula");
    value = 9.0;
    ASSERT_TRUE (&ctx, ml_school_circle_area (-1.0, &value) == ML_ERR_INVALID_ARG && value == 9.0, "negative radius rejected atomically");
    ASSERT_TRUE (&ctx, ml_school_rectangle_area (1e308, 1e308, &value) == ML_ERR_OVERFLOW && value == 9.0,
                 "unrepresentable geometry result rejected atomically");

    return ml_test_summary (&ctx);
}
