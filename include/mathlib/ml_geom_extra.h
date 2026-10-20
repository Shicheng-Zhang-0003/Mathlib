#ifndef MATHLIB_ML_GEOM_EXTRA_H
#define MATHLIB_ML_GEOM_EXTRA_H

#include "ml_compiler.h"
#include "ml_core.h"
#include "ml_types.h"
#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * v12R2 COMPUTATIONAL GEOMETRY & CURVES PACK
 *
 * Minimal enclosing circle, Delaunay triangulation (Bowyer-Watson), Voronoi
 * cells (half-plane clipping), a k-d tree (build / knn / radius), Bezier,
 * B-spline and Catmull-Rom curves, conic classification, Sutherland-Hodgman
 * polygon clipping and 3-D segment-segment closest points.
 * ========================================================================== */

/* ------------------------------------------------- enclosing circle ------ */
ML_API ml_status_t ml_min_enclosing_circle (const double *x, const double *y,
                                            int n, double *cx, double *cy,
                                            double *r);

/* ------------------------------------------------------- Delaunay -------- */
/* Triangulates n points (x,y).  tris receives 3 indices per triangle and must
 * have room for 3*max_tris entries.  Returns the number of triangles, or -1
 * on invalid input. */
ML_API int ml_delaunay_triangulate (const double *x, const double *y, int n,
                                    int *tris, int max_tris);

/* ------------------------------------------------------- Voronoi --------- */
/* Computes the Voronoi cell of point `i` clipped to the axis box
 * [xmin,xmax]x[ymin,ymax].  Writes up to maxv vertices.  Returns the vertex
 * count (0 if the point lies outside the box). */
ML_API int ml_voronoi_cell (int i, const double *x, const double *y, int n,
                            double xmin, double xmax, double ymin, double ymax,
                            double *outx, double *outy, int maxv);
/* --------------------------------------------------------- k-d tree ------ */
typedef struct {
    int point; /* index into the point array; -1 marks an empty node */
    int left;
    int right;
    int axis;
} ml_kdnode_t;
/* Builds a k-d tree over n points of dimension dim.  nodes must hold at least
 * n entries.  Returns the root node index, or -1 on error. */
ML_API int ml_kdtree_build (const double *pts, int n, int dim,
                            ml_kdnode_t *nodes, int max_nodes);

/* k nearest neighbours of query q.  Writes k indices (ascending distance) to
 * out_idx and squared distances to out_dist.  Returns the number found. */
ML_API int ml_kdtree_knn (const double *pts, const ml_kdnode_t *nodes, int root,
                          int dim, const double *q, int k, int *out_idx,
                          double *out_dist);

/* All points within `radius` of q.  Writes indices to out_idx.  Returns the
 * number found (bounded by max_out). */
ML_API int ml_kdtree_radius (const double *pts, const ml_kdnode_t *nodes,
                             int root, int dim, const double *q, double radius,
                             int *out_idx, int max_out);

/* ---------------------------------------------------------- curves ------- */
ML_API void ml_bezier_eval (const double *ctrl, int n, int dim, double t,
                            double *out);
ML_API void ml_bezier_deriv (const double *ctrl, int n, int dim, double t,
                             double *out);
ML_API void ml_bspline_eval (const double *ctrl, int nctrl, int degree,
                             const double *knots, double u, int dim,
                             double *out);
ML_API void ml_catmull_rom (const double *p0, const double *p1,
                            const double *p2, const double *p3, double t,
                            int dim, double *out);
/* ---------------------------------------------------------- conics ------- */
typedef enum {
    ML_CONIC_EMPTY = 0,
    ML_CONIC_ELLIPSE = 1,
    ML_CONIC_CIRCLE = 2,
    ML_CONIC_PARABOLA = 3,
    ML_CONIC_HYPERBOLA = 4,
    ML_CONIC_LINE = 5,
    ML_CONIC_POINT = 6
} ml_conic_type_t;
/* General conic A x^2 + B xy + C y^2 + D x + E y + F = 0. */
ML_API ml_conic_type_t ml_conic_classify (double A, double B, double C,
                                          double D, double E, double F);

/* ------------------------------------------------------- polygons -------- */
/* Sutherland-Hodgman clip of a subject polygon by a convex clip polygon.
 * Returns the output vertex count. */
ML_API int ml_polygon_clip (const double *sx, const double *sy, int n,
                            const double *cx, const double *cy, int m,
                            double *ox, double *oy, int max_out);

/* ------------------------------------------------------------ 3-D -------- */
ML_API void ml_seg_seg_closest (const double *a0, const double *a1,
                                const double *b0, const double *b1,
                                double *pa, double *pb);

#ifdef __cplusplus
}
#endif
#endif /* MATHLIB_ML_GEOM_EXTRA_H */
