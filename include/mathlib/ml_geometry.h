#ifndef MATHLIB_ML_GEOMETRY_H
#define MATHLIB_ML_GEOMETRY_H

#include "ml_compiler.h"
#include "ml_core.h"
#include "ml_types.h"
#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Euclidean / coordinate geometry: triangles, polygons, circles, 3D basics.
 *
 * Fills the high-school pillar (only Heron + SAS existed) and the IMO gap
 * (centers, Euler line, power of a point, exact predicates). Conventions:
 * - Areas/distances return double, NaN on non-finite input or degeneracy.
 * - Multi-output constructors return ml_status_t and NaN-poison outs on
 *   failure (SINGULAR for collinear/degenerate input).
 * - int verdicts (intersection, containment, Euler/nine-point checks) fail
 *   closed to 0 on invalid input. Touching counts as intersecting;
 *   on-edge counts as contained (documented).
 * - Orientation predicates run in long double: exact for exactly
 *   representable coordinates (the contest-integer case).
 * ========================================================================== */

/* Areas */
ML_API double ml_tri_area_coords (double x1, double y1, double x2, double y2,
                                  double x3, double y3);
ML_API double ml_polygon_area (const double *xs, const double *ys, int n);

/* Triangle centers (NaN-poison + SINGULAR when collinear) */
ML_API ml_status_t ml_tri_centroid (double x1, double y1, double x2, double y2,
                                    double x3, double y3, double *ox, double *oy);
ML_API ml_status_t ml_tri_circumcenter (double x1, double y1, double x2, double y2,
                                        double x3, double y3, double *ox, double *oy);
ML_API ml_status_t ml_tri_incenter (double x1, double y1, double x2, double y2,
                                    double x3, double y3, double *ox, double *oy);
ML_API ml_status_t ml_tri_orthocenter (double x1, double y1, double x2, double y2,
                                       double x3, double y3, double *ox, double *oy);
ML_API double ml_tri_circumradius (double a, double b, double c);
ML_API double ml_tri_inradius (double a, double b, double c);

/* Classical verifiers (tol-gated 1/0) */
ML_API int ml_euler_line_check (double x1, double y1, double x2, double y2,
                                double x3, double y3, double tol);
ML_API int ml_nine_point_check (double x1, double y1, double x2, double y2,
                                double x3, double y3, double tol);

/* Predicates */
ML_API double ml_power_of_point (double px, double py, double cx, double cy, double r);
ML_API int ml_seg_intersect (double ax, double ay, double bx, double by,
                             double cx, double cy, double dx, double dy);
ML_API int ml_point_in_polygon (double x, double y, const double *xs,
                                const double *ys, int n);
ML_API int ml_convex_hull (const double *xs, const double *ys, int n,
                           int *out, int cap);

/* Distances: 2D segment here; 3D plane reuses the hardened
 * ml_point_plane_dist in ml_linalg.h (scaled-norm, overflow-safe). */
ML_API double ml_point_seg_dist (double px, double py, double ax, double ay,
                                 double bx, double by);

/* Recreational + solids */
ML_API double ml_clock_angle (double h, double m);
ML_API double ml_sphere_vol (double r);
ML_API double ml_sphere_area (double r);
ML_API double ml_cone_vol (double r, double h);
ML_API double ml_cyl_vol (double r, double h);

#ifdef __cplusplus
}
#endif
#endif /* MATHLIB_ML_GEOMETRY_H */
