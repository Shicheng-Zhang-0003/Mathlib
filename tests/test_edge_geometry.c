/* Geometry pins: centers, radii, Euler/nine-point, predicates, hull.
 * Ground truth: hand-computed exact (3-4-5 triangle, unit square). */
#include "test_harness.h"
#include "ml_geometry.h"
#include "ml_core.h"

int main(void) {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Edge Geometry");

    /* 3-4-5 triangle A=(0,0), B=(4,0), C=(0,3). */
    ASSERT_NEAR(&ctx, ml_tri_area_coords(0, 0, 4, 0, 0, 3), 6.0, 1e-12, "area 3-4-5");
    ASSERT_TRUE(&ctx, ml_isnan(ml_tri_area_coords(0, 0, ml_make_nan(), 0, 0, 3)), "area NaN");
    {
        double ox = 0, oy = 0;
        ASSERT_TRUE(&ctx, ml_tri_centroid(0, 0, 4, 0, 0, 3, &ox, &oy) == ML_SUCCESS, "centroid st");
        ASSERT_NEAR(&ctx, ox, 1.3333333333333333, 1e-12, "centroid x");
        ASSERT_NEAR(&ctx, oy, 1.0, 1e-12, "centroid y");
        ASSERT_TRUE(&ctx, ml_tri_circumcenter(0, 0, 4, 0, 0, 3, &ox, &oy) == ML_SUCCESS, "cc st");
        ASSERT_NEAR(&ctx, ox, 2.0, 1e-12, "cc x");
        ASSERT_NEAR(&ctx, oy, 1.5, 1e-12, "cc y");
        ASSERT_TRUE(&ctx, ml_tri_incenter(0, 0, 4, 0, 0, 3, &ox, &oy) == ML_SUCCESS, "ic st");
        ASSERT_NEAR(&ctx, ox, 1.0, 1e-12, "ic x");
        ASSERT_NEAR(&ctx, oy, 1.0, 1e-12, "ic y");
        ASSERT_TRUE(&ctx, ml_tri_orthocenter(0, 0, 4, 0, 0, 3, &ox, &oy) == ML_SUCCESS, "oc st");
        ASSERT_NEAR(&ctx, ox, 0.0, 1e-12, "oc x");
        ASSERT_NEAR(&ctx, oy, 0.0, 1e-12, "oc y");
        ASSERT_TRUE(&ctx, ml_tri_circumcenter(0, 0, 1, 1, 2, 2, &ox, &oy) == ML_ERR_SINGULAR, "collinear singular");
    }
    ASSERT_NEAR(&ctx, ml_tri_circumradius(3, 4, 5), 2.5, 1e-12, "R 3-4-5");
    ASSERT_NEAR(&ctx, ml_tri_inradius(3, 4, 5), 1.0, 1e-12, "r 3-4-5");
    ASSERT_TRUE(&ctx, ml_isnan(ml_tri_circumradius(1, 2, 3)), "degenerate R NaN");
    ASSERT_TRUE(&ctx, ml_euler_line_check(0, 0, 4, 0, 0, 3, 1e-9) == 1, "euler holds");
    ASSERT_TRUE(&ctx, ml_euler_line_check(0, 0, 1, 0, 0, 1, 1e-9) == 1, "euler right-isosceles");
    ASSERT_TRUE(&ctx, ml_nine_point_check(0, 0, 4, 0, 0, 3, 1e-9) == 1, "nine-point holds");
    ASSERT_TRUE(&ctx, ml_nine_point_check(0, 0, 1, 1, 2, 2, 1e-9) == 0, "nine-point degenerate");

    /* Polygons + predicates. */
    {
        double xs[] = {0, 1, 1, 0}, ys[] = {0, 0, 1, 1};
        ASSERT_NEAR(&ctx, ml_polygon_area(xs, ys, 4), 1.0, 1e-12, "square area");
        ASSERT_TRUE(&ctx, ml_point_in_polygon(0.5, 0.5, xs, ys, 4) == 1, "inside");
        ASSERT_TRUE(&ctx, ml_point_in_polygon(2, 2, xs, ys, 4) == 0, "outside");
        ASSERT_TRUE(&ctx, ml_point_in_polygon(0, 0.5, xs, ys, 4) == 1, "on-edge inside");
    }
    ASSERT_NEAR(&ctx, ml_power_of_point(5, 0, 0, 0, 3), 16.0, 1e-12, "power outside");
    ASSERT_NEAR(&ctx, ml_power_of_point(1, 0, 0, 0, 3), -8.0, 1e-12, "power inside");
    ASSERT_TRUE(&ctx, ml_seg_intersect(0, 0, 2, 2, 0, 2, 2, 0) == 1, "crossing");
    ASSERT_TRUE(&ctx, ml_seg_intersect(0, 0, 1, 0, 0, 1, 1, 1) == 0, "parallel apart");
    ASSERT_TRUE(&ctx, ml_seg_intersect(0, 0, 2, 0, 1, 0, 3, 0) == 1, "collinear overlap");
    ASSERT_TRUE(&ctx, ml_seg_intersect(0, 0, 1, 1, 1, 1, 2, 0) == 1, "endpoint touch");
    {
        double xs[] = {0, 1, 1, 0, 0.5}, ys[] = {0, 0, 1, 1, 0.5};
        int out[8];
        ASSERT_TRUE(&ctx, ml_convex_hull(xs, ys, 5, NULL, 0) == 4, "hull query count");
        ASSERT_TRUE(&ctx, ml_convex_hull(xs, ys, 5, out, 8) == 4, "hull 4 corners");
        ASSERT_TRUE(&ctx, ml_convex_hull(xs, ys, 5, out, 2) == -1, "hull cap short");
    }

    /* Distances, clock, solids. */
    ASSERT_NEAR(&ctx, ml_point_seg_dist(1, 1, 0, 0, 2, 0), 1.0, 1e-12, "pt-seg");
    ASSERT_NEAR(&ctx, ml_point_seg_dist(3, 0, 0, 0, 2, 0), 1.0, 1e-12, "pt-seg clamp");
    ASSERT_NEAR(&ctx, ml_clock_angle(3, 0), 90.0, 1e-12, "clock 3:00");
    ASSERT_NEAR(&ctx, ml_clock_angle(12 - 12, 0), 0.0, 1e-12, "clock 0:00");
    ASSERT_TRUE(&ctx, ml_isnan(ml_clock_angle(12, 0)), "clock invalid");
    ASSERT_NEAR(&ctx, ml_sphere_vol(1), 4.1887902047863905, 1e-12, "sphere vol");
    ASSERT_NEAR(&ctx, ml_sphere_area(1), 12.566370614359172, 1e-12, "sphere area");
    ASSERT_NEAR(&ctx, ml_cone_vol(3, 2), 18.849555921538759, 1e-12, "cone vol");
    ASSERT_NEAR(&ctx, ml_cyl_vol(2, 3), 37.699111843077517, 1e-12, "cyl vol");
    ASSERT_TRUE(&ctx, ml_isnan(ml_sphere_vol(-1)), "solid negative NaN");

    return ml_test_summary(&ctx);
}
