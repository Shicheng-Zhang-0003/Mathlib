#include "ml_compiler.h"
#include "ml_geometry.h"
#include "ml_core.h"
#include "ml_trig.h"
#include <stddef.h>
#include <stdlib.h>
/* Euclidean / coordinate geometry. Predicates in long double (exact for
 * exactly-representable coordinates); constructors NaN-poison on failure. */

static int mlg_finite2 (double x, double y) {
    return ml_isfinite (x) && ml_isfinite (y);
}
static long double mlg_cross2 (long double ax, long double ay,
                               long double bx, long double by) {
    return ax * by - ay * bx;
}
static long double mlg_area2_ld (double x1, double y1, double x2, double y2,
                                 double x3, double y3) {
    return mlg_cross2 ((long double) x2 - (long double) x1,
                       (long double) y2 - (long double) y1,
                       (long double) x3 - (long double) x1,
                       (long double) y3 - (long double) y1);
}
ML_API double ml_tri_area_coords (double x1, double y1, double x2, double y2,
                                  double x3, double y3) {
    long double a2;
    double r;
    if (!mlg_finite2 (x1, y1) || !mlg_finite2 (x2, y2) || !mlg_finite2 (x3, y3))
        return ml_make_nan ();
    a2 = mlg_area2_ld (x1, y1, x2, y2, x3, y3);
    if (a2 < 0.0L) a2 = -a2;
    r = (double) (0.5L * a2);
    return ml_isfinite (r) ? r : ml_make_nan ();
}
ML_API double ml_polygon_area (const double *xs, const double *ys, int n) {
    long double s = 0.0L;
    int i;
    double r;
    if (!xs || !ys || n < 3 || n > 1000000) return ml_make_nan ();
    for (i = 0; i < n; i++) {
        int j = (i + 1) % n;
        if (!mlg_finite2 (xs [i], ys [i])) return ml_make_nan ();
        s += (long double) xs [i] * (long double) ys [j] - (long double) xs [j] * (long double) ys [i];
    }
    if (s < 0.0L) s = -s;
    r = (double) (0.5L * s);
    return ml_isfinite (r) ? r : ml_make_nan ();
}
ML_API ml_status_t ml_tri_centroid (double x1, double y1, double x2, double y2,
                                    double x3, double y3, double *ox, double *oy) {
    if (ML_UNLIKELY (!ox || !oy)) return ML_ERR_INVALID_ARG;
    if (!mlg_finite2 (x1, y1) || !mlg_finite2 (x2, y2) || !mlg_finite2 (x3, y3)) {
        *ox = ml_make_nan ();
        *oy = ml_make_nan ();
        return ML_ERR_NAN_INPUT;
    }
    if (mlg_area2_ld (x1, y1, x2, y2, x3, y3) == 0.0L) {
        *ox = ml_make_nan ();
        *oy = ml_make_nan ();
        return ML_ERR_SINGULAR;
    }
    *ox = (x1 + x2 + x3) / 3.0;
    *oy = (y1 + y2 + y3) / 3.0;
    return ML_SUCCESS;
}
ML_API ml_status_t ml_tri_circumcenter (double x1, double y1, double x2, double y2,
                                        double x3, double y3, double *ox, double *oy) {
    long double d;
    long double s1;
    long double s2;
    long double s3;
    long double ux;
    long double uy;
    long double sc;
    if (ML_UNLIKELY (!ox || !oy)) return ML_ERR_INVALID_ARG;
    if (!mlg_finite2 (x1, y1) || !mlg_finite2 (x2, y2) || !mlg_finite2 (x3, y3)) {
        *ox = ml_make_nan ();
        *oy = ml_make_nan ();
        return ML_ERR_NAN_INPUT;
    }
    d = mlg_area2_ld (x1, y1, x2, y2, x3, y3);
    sc = ml_fabs (x1) + ml_fabs (y1) + ml_fabs (x2) + ml_fabs (y2) + ml_fabs (x3) + ml_fabs (y3);
    if (!(ml_fabs ((double) d) > 1e-15 * (1.0 + sc) * (1.0 + sc))) {
        *ox = ml_make_nan ();
        *oy = ml_make_nan ();
        return ML_ERR_SINGULAR;
    }
    s1 = (long double) x1 * (long double) x1 + (long double) y1 * (long double) y1;
    s2 = (long double) x2 * (long double) x2 + (long double) y2 * (long double) y2;
    s3 = (long double) x3 * (long double) x3 + (long double) y3 * (long double) y3;
    ux = (s1 * ((long double) y2 - (long double) y3) + s2 * ((long double) y3 - (long double) y1) + s3 * ((long double) y1 - (long double) y2)) / (2.0L * d);
    uy = (s1 * ((long double) x3 - (long double) x2) + s2 * ((long double) x1 - (long double) x3) + s3 * ((long double) x2 - (long double) x1)) / (2.0L * d);
    *ox = (double) ux;
    *oy = (double) uy;
    if (!ml_isfinite (*ox) || !ml_isfinite (*oy)) {
        *ox = ml_make_nan ();
        *oy = ml_make_nan ();
        return ML_ERR_OVERFLOW;
    }
    return ML_SUCCESS;
}
ML_API ml_status_t ml_tri_incenter (double x1, double y1, double x2, double y2,
                                    double x3, double y3, double *ox, double *oy) {
    /* Weighted by side lengths: (a*A + b*B + c*C)/(a+b+c), a = |BC|. */
    double a;
    double b;
    double c;
    double per;
    if (ML_UNLIKELY (!ox || !oy)) return ML_ERR_INVALID_ARG;
    if (!mlg_finite2 (x1, y1) || !mlg_finite2 (x2, y2) || !mlg_finite2 (x3, y3)) {
        *ox = ml_make_nan ();
        *oy = ml_make_nan ();
        return ML_ERR_NAN_INPUT;
    }
    if (mlg_area2_ld (x1, y1, x2, y2, x3, y3) == 0.0L) {
        *ox = ml_make_nan ();
        *oy = ml_make_nan ();
        return ML_ERR_SINGULAR;
    }
    a = ml_hypot (x2 - x3, y2 - y3);
    b = ml_hypot (x1 - x3, y1 - y3);
    c = ml_hypot (x1 - x2, y1 - y2);
    per = a + b + c;
    if (!(per > 0.0) || !ml_isfinite (per)) {
        *ox = ml_make_nan ();
        *oy = ml_make_nan ();
        return ML_ERR_SINGULAR;
    }
    *ox = (a * x1 + b * x2 + c * x3) / per;
    *oy = (a * y1 + b * y2 + c * y3) / per;
    return ML_SUCCESS;
}
ML_API ml_status_t ml_tri_orthocenter (double x1, double y1, double x2, double y2,
                                       double x3, double y3, double *ox, double *oy) {
    /* Euler vector relation: H = A + B + C - 2O. */
    double cx;
    double cy;
    ml_status_t st;
    if (ML_UNLIKELY (!ox || !oy)) return ML_ERR_INVALID_ARG;
    st = ml_tri_circumcenter (x1, y1, x2, y2, x3, y3, &cx, &cy);
    if (st != ML_SUCCESS) {
        *ox = ml_make_nan ();
        *oy = ml_make_nan ();
        return st;
    }
    *ox = x1 + x2 + x3 - 2.0 * cx;
    *oy = y1 + y2 + y3 - 2.0 * cy;
    if (!ml_isfinite (*ox) || !ml_isfinite (*oy)) {
        *ox = ml_make_nan ();
        *oy = ml_make_nan ();
        return ML_ERR_OVERFLOW;
    }
    return ML_SUCCESS;
}
ML_API double ml_tri_circumradius (double a, double b, double c) {
    /* R = abc/4K via staged Heron (reuses the validated kernel). */
    double k;
    long double r;
    double o;
    if (ml_isnan (a) || ml_isnan (b) || ml_isnan (c)) return ml_make_nan ();
    if (!(a > 0.0) || !(b > 0.0) || !(c > 0.0)) return ml_make_nan ();
    if (!ml_isfinite (a) || !ml_isfinite (b) || !ml_isfinite (c)) return ml_make_nan ();
    k = ml_heron (a, b, c);
    if (!(k > 0.0) || !ml_isfinite (k)) return ml_make_nan ();
    r = (long double) a * (long double) b * (long double) c / (4.0L * (long double) k);
    o = (double) r;
    return ml_isfinite (o) ? o : ml_make_nan ();
}
ML_API double ml_tri_inradius (double a, double b, double c) {
    /* r = K/s. */
    double k;
    double s;
    double o;
    if (ml_isnan (a) || ml_isnan (b) || ml_isnan (c)) return ml_make_nan ();
    if (!(a > 0.0) || !(b > 0.0) || !(c > 0.0)) return ml_make_nan ();
    if (!ml_isfinite (a) || !ml_isfinite (b) || !ml_isfinite (c)) return ml_make_nan ();
    k = ml_heron (a, b, c);
    if (!(k > 0.0) || !ml_isfinite (k)) return ml_make_nan ();
    s = 0.5 * (a + b + c);
    if (!(s > 0.0) || !ml_isfinite (s)) return ml_make_nan ();
    o = k / s;
    return ml_isfinite (o) ? o : ml_make_nan ();
}
ML_API int ml_euler_line_check (double x1, double y1, double x2, double y2,
                                double x3, double y3, double tol) {
    /* O, G, H collinear with OG:GH = 1:2. Scale-aware tol (default 1e-9). */
    double ox, oy, gx, gy, hx, hy;
    long double ogx, ogy, ghx, ghy;
    long double og;
    long double gh;
    long double cross;
    long double ratio;
    if (!(tol > 0.0) || !ml_isfinite (tol)) tol = 1e-9;
    if (ml_tri_circumcenter (x1, y1, x2, y2, x3, y3, &ox, &oy) != ML_SUCCESS) return 0;
    if (ml_tri_centroid (x1, y1, x2, y2, x3, y3, &gx, &gy) != ML_SUCCESS) return 0;
    if (ml_tri_orthocenter (x1, y1, x2, y2, x3, y3, &hx, &hy) != ML_SUCCESS) return 0;
    ogx = (long double) gx - (long double) ox;
    ogy = (long double) gy - (long double) oy;
    ghx = (long double) hx - (long double) gx;
    ghy = (long double) hy - (long double) gy;
    og = __builtin_sqrtl (ogx * ogx + ogy * ogy);
    gh = __builtin_sqrtl (ghx * ghx + ghy * ghy);
    if (!(og > 0.0L) || !(gh > 0.0L)) {
        /* Equilateral: O = G = H coincides; line degenerates but holds. */
        return (og == 0.0L && gh == 0.0L) ? 1 : 0;
    }
    cross = ogx * ghy - ogy * ghx;
    if (ml_fabs ((double) (cross / (og * gh))) > tol) return 0;
    ratio = og * 2.0L / gh;
    if (ml_fabs ((double) (ratio - 1.0L)) > tol * 10.0L) return 0;
    return 1;
}
ML_API int ml_nine_point_check (double x1, double y1, double x2, double y2,
                                double x3, double y3, double tol) {
    /* N = midpoint(O,H); dist(N, each side midpoint) = R/2 within tol. */
    double ox, oy, hx, hy;
    double nx;
    double ny;
    double mx [3];
    double my [3];
    double a;
    double b;
    double c;
    double r;
    int i;
    if (!(tol > 0.0) || !ml_isfinite (tol)) tol = 1e-9;
    if (ml_tri_circumcenter (x1, y1, x2, y2, x3, y3, &ox, &oy) != ML_SUCCESS) return 0;
    if (ml_tri_orthocenter (x1, y1, x2, y2, x3, y3, &hx, &hy) != ML_SUCCESS) return 0;
    nx = 0.5 * (ox + hx);
    ny = 0.5 * (oy + hy);
    if (!ml_isfinite (nx) || !ml_isfinite (ny)) return 0;
    mx [0] = 0.5 * (x1 + x2);
    my [0] = 0.5 * (y1 + y2);
    mx [1] = 0.5 * (x2 + x3);
    my [1] = 0.5 * (y2 + y3);
    mx [2] = 0.5 * (x3 + x1);
    my [2] = 0.5 * (y3 + y1);
    a = ml_hypot (x2 - x3, y2 - y3);
    b = ml_hypot (x1 - x3, y1 - y3);
    c = ml_hypot (x1 - x2, y1 - y2);
    r = ml_tri_circumradius (a, b, c);
    if (!(r > 0.0) || !ml_isfinite (r)) return 0;
    for (i = 0; i < 3; i++) {
        double d = ml_hypot (mx [i] - nx, my [i] - ny);
        if (ml_fabs (d - 0.5 * r) > tol * (1.0 + r)) return 0;
    }
    return 1;
}
ML_API double ml_power_of_point (double px, double py, double cx, double cy, double r) {
    long double dx;
    long double dy;
    double o;
    if (!mlg_finite2 (px, py) || !mlg_finite2 (cx, cy)) return ml_make_nan ();
    if (ml_isnan (r) || !ml_isfinite (r) || r < 0.0) return ml_make_nan ();
    dx = (long double) px - (long double) cx;
    dy = (long double) py - (long double) cy;
    o = (double) (dx * dx + dy * dy - (long double) r * (long double) r);
    return ml_isfinite (o) ? o : ml_make_nan ();
}
static int mlg_orient_sign (double ax, double ay, double bx, double by,
                            double cx, double cy) {
    long double v = mlg_cross2 ((long double) bx - (long double) ax,
                                (long double) by - (long double) ay,
                                (long double) cx - (long double) ax,
                                (long double) cy - (long double) ay);
    if (v > 0.0L) return 1;
    if (v < 0.0L) return -1;
    return 0;
}
static int mlg_on_seg (double ax, double ay, double bx, double by,
                       double cx, double cy) {
    /* B on segment AC (collinear caller-side). Long-double exact. */
    long double minx = (long double) ax < (long double) cx ? (long double) ax : (long double) cx;
    long double maxx = (long double) ax > (long double) cx ? (long double) ax : (long double) cx;
    long double miny = (long double) ay < (long double) cy ? (long double) ay : (long double) cy;
    long double maxy = (long double) ay > (long double) cy ? (long double) ay : (long double) cy;
    return ((long double) bx >= minx && (long double) bx <= maxx &&
            (long double) by >= miny && (long double) by <= maxy)
               ? 1
               : 0;
}
ML_API int ml_seg_intersect (double ax, double ay, double bx, double by,
                             double cx, double cy, double dx, double dy) {
    int o1;
    int o2;
    int o3;
    int o4;
    if (!mlg_finite2 (ax, ay) || !mlg_finite2 (bx, by)) return 0;
    if (!mlg_finite2 (cx, cy) || !mlg_finite2 (dx, dy)) return 0;
    o1 = mlg_orient_sign (ax, ay, bx, by, cx, cy);
    o2 = mlg_orient_sign (ax, ay, bx, by, dx, dy);
    o3 = mlg_orient_sign (cx, cy, dx, dy, ax, ay);
    o4 = mlg_orient_sign (cx, cy, dx, dy, bx, by);
    if (o1 != o2 && o3 != o4) return 1;
    if (o1 == 0 && mlg_on_seg (ax, ay, cx, cy, bx, by)) return 1;
    if (o2 == 0 && mlg_on_seg (ax, ay, dx, dy, bx, by)) return 1;
    if (o3 == 0 && mlg_on_seg (cx, cy, ax, ay, dx, dy)) return 1;
    if (o4 == 0 && mlg_on_seg (cx, cy, bx, by, dx, dy)) return 1;
    return 0;
}
ML_API int ml_point_in_polygon (double x, double y, const double *xs,
                                const double *ys, int n) {
    int inside = 0;
    int i;
    if (ml_isnan (x) || ml_isnan (y) || !ml_isfinite (x) || !ml_isfinite (y)) return 0;
    if (!xs || !ys || n < 3 || n > 1000000) return 0;
    for (i = 0; i < n; i++) {
        int j = (i + 1) % n;
        double xi = xs [i], yi = ys [i], xj = xs [j], yj = ys [j];
        long double cross;
        if (!mlg_finite2 (xi, yi)) return 0;
        /* On-edge (exact LD) counts as contained. */
        cross = mlg_cross2 ((long double) xj - (long double) xi,
                            (long double) yj - (long double) yi,
                            (long double) x - (long double) xi,
                            (long double) y - (long double) yi);
        if (cross == 0.0L && mlg_on_seg (xi, yi, x, y, xj, yj)) return 1;
        if (((yi > y) != (yj > y))) {
            long double xint = (long double) xi + ((long double) y - (long double) yi) /
                                                      ((long double) yj - (long double) yi) * ((long double) xj - (long double) xi);
            if ((long double) x < xint) inside = !inside;
        }
    }
    return inside;
}
typedef struct {
    double x;
    double y;
    int idx;
} mlg_pt_t;
static int mlg_pt_cmp (const void *pa, const void *pb) {
    const mlg_pt_t *a = (const mlg_pt_t *) pa;
    const mlg_pt_t *b = (const mlg_pt_t *) pb;
    if (a->x < b->x) return -1;
    if (a->x > b->x) return 1;
    if (a->y < b->y) return -1;
    if (a->y > b->y) return 1;
    return 0;
}
ML_API int ml_convex_hull (const double *xs, const double *ys, int n,
                           int *out, int cap) {
    /* Andrew monotone chain. out == NULL queries the count. Collinear
     * points on edges are dropped (extreme endpoints kept). */
    mlg_pt_t *pts = NULL;
    int *hull = NULL;
    int m = 0;
    int i;
    int k = 0;
    if (!xs || !ys || n < 1 || n > 100000) return -1;
    if ((size_t) n > (size_t) -1 / sizeof (mlg_pt_t)) return -1;
    pts = (mlg_pt_t *) malloc ((size_t) n * sizeof (mlg_pt_t));
    if (!pts) return -1;
    for (i = 0; i < n; i++) {
        if (!mlg_finite2 (xs [i], ys [i])) {
            free (pts);
            return -1;
        }
        pts [i].x = xs [i];
        pts [i].y = ys [i];
        pts [i].idx = i;
    }
    qsort (pts, (size_t) n, sizeof (mlg_pt_t), mlg_pt_cmp);
    /* Dedup. */
    m = 0;
    for (i = 0; i < n; i++) {
        if (m > 0 && pts [i].x == pts [m - 1].x && pts [i].y == pts [m - 1].y) continue;
        pts [m++] = pts [i];
    }
    if (m == 1) {
        if (out && cap >= 1) out [0] = pts [0].idx;
        free (pts);
        return (out && cap < 1) ? -1 : 1;
    }
    if ((size_t) (2 * m) > (size_t) -1 / sizeof (int)) {
        free (pts);
        return -1;
    }
    hull = (int *) malloc ((size_t) (2 * m) * sizeof (int));
    if (!hull) {
        free (pts);
        return -1;
    }
    for (i = 0; i < m; i++) {
        while (k >= 2 && mlg_orient_sign (pts [hull [k - 2]].x, pts [hull [k - 2]].y,
                                          pts [hull [k - 1]].x, pts [hull [k - 1]].y,
                                          pts [i].x, pts [i].y) <= 0) k--;
        hull [k++] = i;
    }
    {
        int t = k + 1;
        for (i = m - 2; i >= 0; i--) {
            while (k >= t && mlg_orient_sign (pts [hull [k - 2]].x, pts [hull [k - 2]].y,
                                              pts [hull [k - 1]].x, pts [hull [k - 1]].y,
                                              pts [i].x, pts [i].y) <= 0) k--;
            hull [k++] = i;
        }
    }
    k--;
    if (out == NULL) {
        free (hull);
        free (pts);
        return k;
    }
    if (cap < k) {
        free (hull);
        free (pts);
        return -1;
    }
    for (i = 0; i < k; i++) out [i] = pts [hull [i]].idx;
    free (hull);
    free (pts);
    return k;
}
ML_API double ml_point_seg_dist (double px, double py, double ax, double ay,
                                 double bx, double by) {
    long double abx;
    long double aby;
    long double denom;
    long double t;
    long double cx;
    long double cy;
    double o;
    if (!mlg_finite2 (px, py) || !mlg_finite2 (ax, ay) || !mlg_finite2 (bx, by))
        return ml_make_nan ();
    abx = (long double) bx - (long double) ax;
    aby = (long double) by - (long double) ay;
    denom = abx * abx + aby * aby;
    if (denom == 0.0L) {
        o = ml_hypot (px - ax, py - ay);
        return ml_isfinite (o) ? o : ml_make_nan ();
    }
    t = ((long double) px - (long double) ax) * abx + ((long double) py - (long double) ay) * aby;
    t /= denom;
    if (t < 0.0L) t = 0.0L;
    if (t > 1.0L) t = 1.0L;
    cx = (long double) ax + t * abx - (long double) px;
    cy = (long double) ay + t * aby - (long double) py;
    o = (double) __builtin_sqrtl (cx * cx + cy * cy);
    return ml_isfinite (o) ? o : ml_make_nan ();
}
ML_API double ml_clock_angle (double h, double m) {
    double a;
    if (ml_isnan (h) || ml_isnan (m) || !ml_isfinite (h) || !ml_isfinite (m)) return ml_make_nan ();
    if (h < 0.0 || h >= 12.0 || m < 0.0 || m >= 60.0) return ml_make_nan ();
    a = ml_fabs (30.0 * h - 5.5 * m);
    if (a > 180.0) a = 360.0 - a;
    return a;
}
ML_API double ml_sphere_vol (double r) {
    double o;
    if (ml_isnan (r) || !ml_isfinite (r) || r < 0.0) return ml_make_nan ();
    o = (4.0 / 3.0) * ML_PI * r * r * r;
    return ml_isfinite (o) ? o : ml_make_nan ();
}
ML_API double ml_sphere_area (double r) {
    double o;
    if (ml_isnan (r) || !ml_isfinite (r) || r < 0.0) return ml_make_nan ();
    o = 4.0 * ML_PI * r * r;
    return ml_isfinite (o) ? o : ml_make_nan ();
}
ML_API double ml_cone_vol (double r, double h) {
    double o;
    if (ml_isnan (r) || ml_isnan (h)) return ml_make_nan ();
    if (!ml_isfinite (r) || !ml_isfinite (h) || r < 0.0 || h < 0.0) return ml_make_nan ();
    o = ML_PI * r * r * h / 3.0;
    return ml_isfinite (o) ? o : ml_make_nan ();
}
ML_API double ml_cyl_vol (double r, double h) {
    double o;
    if (ml_isnan (r) || ml_isnan (h)) return ml_make_nan ();
    if (!ml_isfinite (r) || !ml_isfinite (h) || r < 0.0 || h < 0.0) return ml_make_nan ();
    o = ML_PI * r * r * h;
    return ml_isfinite (o) ? o : ml_make_nan ();
}
