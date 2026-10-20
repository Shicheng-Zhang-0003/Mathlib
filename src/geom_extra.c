#include "ml_compiler.h"
#include "ml_geom_extra.h"

/* ============================================================================
 * v12R2 computational geometry & curves pack.
 * ========================================================================== */

#define GE_EPS 1e-12
static double ge_dist2 (const double *a, const double *b, int dim) {
    double s = 0.0;
    int i;
    for (i = 0; i < dim; i++) {
        double d = a [i] - b [i];
        s += d * d;
    }
    return s;
}
/* --------------------------------------------- minimal enclosing circle -- */

static int ge_in_circle (const double *x, const double *y, int i, double cx,
                         double cy, double r) {
    double dx = x [i] - cx, dy = y [i] - cy;
    return dx * dx + dy * dy <= r * r + 1e-12;
}
static void ge_circumcircle (const double *x, const double *y, int a, int b,
                             int c, double *cx, double *cy, double *r) {
    double ax = x [a], ay = y [a], bx = x [b], by = y [b], cxx = x [c], cyy = y [c];
    double d = 2.0 * (ax * (by - cyy) + bx * (cyy - ay) + cxx * (ay - by));
    if (ml_fabs (d) < 1e-15) {
        *cx = (ax + bx) / 2.0;
        *cy = (ay + by) / 2.0;
    } else {
        double a2 = ax * ax + ay * ay, b2 = bx * bx + by * by, c2 = cxx * cxx + cyy * cyy;
        *cx = (a2 * (by - cyy) + b2 * (cyy - ay) + c2 * (ay - by)) / d;
        *cy = (a2 * (cxx - bx) + b2 * (ax - cxx) + c2 * (bx - ax)) / d;
    }
    {
        double dx = ax - *cx, dy = ay - *cy;
        *r = ml_sqrt (dx * dx + dy * dy);
    }
}
ML_API ml_status_t ml_min_enclosing_circle (const double *x, const double *y,
                                            int n, double *cx, double *cy,
                                            double *r) {
    double ccx = 0.0, ccy = 0.0, rr = -1.0;
    int i, j, k;
    if (ML_UNLIKELY (!x || !y || !cx || !cy || !r || n <= 0))
        return ML_ERR_INVALID_ARG;
    for (i = 0; i < n; i++) {
        if (ge_in_circle (x, y, i, ccx, ccy, rr)) continue;
        ccx = x [i];
        ccy = y [i];
        rr = 0.0;
        for (j = 0; j < i; j++) {
            if (ge_in_circle (x, y, j, ccx, ccy, rr)) continue;
            ccx = (x [i] + x [j]) / 2.0;
            ccy = (y [i] + y [j]) / 2.0;
            rr = ml_sqrt ((x [i] - ccx) * (x [i] - ccx) + (y [i] - ccy) * (y [i] - ccy));
            for (k = 0; k < j; k++) {
                if (ge_in_circle (x, y, k, ccx, ccy, rr)) continue;
                ge_circumcircle (x, y, i, j, k, &ccx, &ccy, &rr);
            }
        }
    }
    *cx = ccx;
    *cy = ccy;
    *r = rr;
    return ML_SUCCESS;
}
/* ------------------------------------------------------- Delaunay -------- */

#define GE_MAXP 512
#define GE_MAXT 2048
static int ge_orient (const double *px, const double *py, int a, int b, int c) {
    double v = (px [b] - px [a]) * (py [c] - py [a]) - (py [b] - py [a]) * (px [c] - px [a]);
    return v > 0.0 ? 1 : (v < 0.0 ? -1 : 0);
}
static int ge_in_circum (const double *px, const double *py, int a, int b, int c,
                         int p) {
    double ax = px [a] - px [p], ay = py [a] - py [p];
    double bx = px [b] - px [p], by = py [b] - py [p];
    double cx = px [c] - px [p], cy = py [c] - py [p];
    double det = (ax * ax + ay * ay) * (bx * cy - by * cx) - (bx * bx + by * by) * (ax * cy - ay * cx) + (cx * cx + cy * cy) * (ax * by - ay * bx);
    if (ge_orient (px, py, a, b, c) < 0) det = -det;
    return det > 1e-12;
}
ML_API int ml_delaunay_triangulate (const double *x, const double *y, int n,
                                    int *tris, int max_tris) {
    double px [GE_MAXP + 3], py [GE_MAXP + 3];
    int tri [GE_MAXT * 3], ntri = 0;
    int i, t;
    if (ML_UNLIKELY (!x || !y || !tris || n < 3 || n > GE_MAXP))
        return -1;
    {
        double xmin = x [0], xmax = x [0], ymin = y [0], ymax = y [0];
        double dx, dy, dmax, mx, my;
        for (i = 1; i < n; i++) {
            if (x [i] < xmin) xmin = x [i];
            if (x [i] > xmax) xmax = x [i];
            if (y [i] < ymin) ymin = y [i];
            if (y [i] > ymax) ymax = y [i];
        }
        mx = (xmin + xmax) / 2.0;
        my = (ymin + ymax) / 2.0;
        dx = xmax - xmin;
        dy = ymax - ymin;
        dmax = (dx > dy ? dx : dy);
        if (dmax <= 0.0) dmax = 1.0;
        for (i = 0; i < n; i++) {
            px [i] = x [i];
            py [i] = y [i];
        }
        px [n] = mx - 20.0 * dmax;
        py [n] = my - dmax;
        px [n + 1] = mx;
        py [n + 1] = my + 20.0 * dmax;
        px [n + 2] = mx + 20.0 * dmax;
        py [n + 2] = my - dmax;
    }
    tri [0] = n;
    tri [1] = n + 1;
    tri [2] = n + 2;
    ntri = 1;
    for (i = 0; i < n; i++) {
        int bad [GE_MAXT], nbad = 0;
        int edges [GE_MAXT * 2], ne = 0;
        unsigned char isbad [GE_MAXT];
        int q, r, e, w;
        for (t = 0; t < ntri; t++)
            if (ge_in_circum (px, py, tri [t * 3], tri [t * 3 + 1], tri [t * 3 + 2], i))
                bad [nbad++] = t;
        for (t = 0; t < ntri; t++) isbad [t] = 0;
        for (t = 0; t < nbad; t++) isbad [bad [t]] = 1;
        for (q = 0; q < nbad; q++) {
            int tv = bad [q];
            int vv0 = tri [tv * 3], vv1 = tri [tv * 3 + 1], vv2 = tri [tv * 3 + 2];
            int ve [3][2];
            ve [0][0] = vv0;
            ve [0][1] = vv1;
            ve [1][0] = vv1;
            ve [1][1] = vv2;
            ve [2][0] = vv2;
            ve [2][1] = vv0;
            for (e = 0; e < 3; e++) {
                int u1 = ve [e][0], u2 = ve [e][1], shared = 0;
                for (r = 0; r < nbad && !shared; r++) {
                    int rv, rw [3];
                    if (r == q) continue;
                    rv = bad [r];
                    rw [0] = tri [rv * 3];
                    rw [1] = tri [rv * 3 + 1];
                    rw [2] = tri [rv * 3 + 2];
                    for (w = 0; w < 3; w++)
                        if ((rw [w] == u1 && rw [(w + 1) % 3] == u2) ||
                            (rw [w] == u2 && rw [(w + 1) % 3] == u1)) {
                            shared = 1;
                            break;
                        }
                }
                if (!shared) {
                    edges [ne * 2] = u1;
                    edges [ne * 2 + 1] = u2;
                    ne++;
                }
            }
        }
        {
            int wi = 0;
            for (t = 0; t < ntri; t++)
                if (!isbad [t]) {
                    tri [wi * 3] = tri [t * 3];
                    tri [wi * 3 + 1] = tri [t * 3 + 1];
                    tri [wi * 3 + 2] = tri [t * 3 + 2];
                    wi++;
                }
            ntri = wi;
        }
        for (t = 0; t < ne; t++) {
            int u1 = edges [t * 2], u2 = edges [t * 2 + 1];
            if (ge_orient (px, py, u1, u2, i) < 0) {
                int s = u1;
                u1 = u2;
                u2 = s;
            }
            if (ntri >= GE_MAXT) return -1;
            tri [ntri * 3] = u1;
            tri [ntri * 3 + 1] = u2;
            tri [ntri * 3 + 2] = i;
            ntri++;
        }
    }
    {
        int out = 0;
        for (t = 0; t < ntri; t++) {
            int a = tri [t * 3], b = tri [t * 3 + 1], c = tri [t * 3 + 2];
            if (a < n && b < n && c < n && out < max_tris) {
                tris [out * 3] = a;
                tris [out * 3 + 1] = b;
                tris [out * 3 + 2] = c;
                out++;
            }
        }
        return out;
    }
}
/* ------------------------------------------------------- Voronoi --------- */

static int ge_clip_half (const double *px, const double *py, int np,
                         double nx, double ny, double c, double *ox, double *oy) {
    int i, out = 0;
    for (i = 0; i < np; i++) {
        int k = (i + 1) % np;
        double si = nx * px [i] + ny * py [i] - c;
        double sk = nx * px [k] + ny * py [k] - c;
        if (si <= 0.0) {
            ox [out] = px [i];
            oy [out] = py [i];
            out++;
        }
        if ((si < 0.0 && sk > 0.0) || (si > 0.0 && sk < 0.0)) {
            double tt = si / (si - sk);
            ox [out] = px [i] + tt * (px [k] - px [i]);
            oy [out] = py [i] + tt * (py [k] - py [i]);
            out++;
        }
    }
    return out;
}
ML_API int ml_voronoi_cell (int i, const double *x, const double *y, int n,
                            double xmin, double xmax, double ymin, double ymax,
                            double *outx, double *outy, int maxv) {
    double px [256], py [256], qx [256], qy [256];
    int np = 4, j;
    if (ML_UNLIKELY (!x || !y || !outx || !outy || i < 0 || i >= n))
        return 0;
    if (n > 256 || maxv < 4) return 0;
    px [0] = xmin;
    py [0] = ymin;
    px [1] = xmax;
    py [1] = ymin;
    px [2] = xmax;
    py [2] = ymax;
    px [3] = xmin;
    py [3] = ymax;
    for (j = 0; j < n; j++) {
        double nx, ny, c, len2;
        if (j == i) continue;
        nx = x [j] - x [i];
        ny = y [j] - y [i];
        len2 = nx * nx + ny * ny;
        if (len2 <= 0.0) continue;
        c = nx * (x [i] + x [j]) / 2.0 + ny * (y [i] + y [j]) / 2.0;
        np = ge_clip_half (px, py, np, nx, ny, c, qx, qy);
        if (np <= 0) return 0;
        {
            int t;
            for (t = 0; t < np; t++) {
                px [t] = qx [t];
                py [t] = qy [t];
            }
        }
    }
    if (np > maxv) np = maxv;
    for (j = 0; j < np; j++) {
        outx [j] = px [j];
        outy [j] = py [j];
    }
    return np;
}
/* --------------------------------------------------------- k-d tree ------ */

static void ge_swap_int (int *a, int *b) {
    int t = *a;
    *a = *b;
    *b = t;
}
static int ge_partition (const double *pts, int dim, int *idx, int lo, int hi,
                         int axis) {
    double piv = pts [(size_t) idx [hi] * (size_t) dim + (size_t) axis];
    int i = lo, j;
    for (j = lo; j < hi; j++)
        if (pts [(size_t) idx [j] * (size_t) dim + (size_t) axis] < piv)
            ge_swap_int (&idx [i++], &idx [j]);
    ge_swap_int (&idx [i], &idx [hi]);
    return i;
}
static void ge_select (const double *pts, int dim, int *idx, int lo, int hi,
                       int k, int axis) {
    while (lo < hi) {
        int p = ge_partition (pts, dim, idx, lo, hi, axis);
        if (p == k) return;
        if (k < p) hi = p - 1;
        else lo = p + 1;
    }
}
static int ge_build (const double *pts, int n, int dim, ml_kdnode_t *nodes,
                     int max_nodes, int *idx, int lo, int hi, int *counter) {
    int axis = 0, mid, node, i, d;
    double best = -1.0;
    if (lo >= hi) return -1;
    for (d = 0; d < dim; d++) {
        double mn = pts [(size_t) idx [lo] * (size_t) dim + (size_t) d];
        double mx = mn;
        for (i = lo + 1; i < hi; i++) {
            double v = pts [(size_t) idx [i] * (size_t) dim + (size_t) d];
            if (v < mn) mn = v;
            if (v > mx) mx = v;
        }
        if (mx - mn > best) {
            best = mx - mn;
            axis = d;
        }
    }
    mid = lo + (hi - lo) / 2;
    ge_select (pts, dim, idx, lo, hi - 1, mid, axis);
    node = (*counter)++;
    if (node >= max_nodes) return -1;
    nodes [node].point = idx [mid];
    nodes [node].axis = axis;
    nodes [node].left = ge_build (pts, n, dim, nodes, max_nodes, idx, lo, mid, counter);
    nodes [node].right = ge_build (pts, n, dim, nodes, max_nodes, idx, mid + 1, hi, counter);
    return node;
}
ML_API int ml_kdtree_build (const double *pts, int n, int dim,
                            ml_kdnode_t *nodes, int max_nodes) {
    int idx [GE_MAXP], i, counter = 0;
    if (ML_UNLIKELY (!pts || !nodes || n <= 0 || dim <= 0 || n > GE_MAXP))
        return -1;
    for (i = 0; i < n; i++) idx [i] = i;
    return ge_build (pts, n, dim, nodes, max_nodes, idx, 0, n, &counter);
}
typedef struct {
    const double *pts;
    int dim;
    double q [16];
    int k;
    int idx [64];
    double dist [64];
    int count;
} ge_knn_t;
static void ge_knn_insert (ge_knn_t *s, int point, double d2) {
    int i;
    if (s->count < s->k) {
        i = s->count++;
    } else if (d2 < s->dist [s->count - 1]) {
        i = s->count - 1;
    } else {
        return;
    }
    s->idx [i] = point;
    s->dist [i] = d2;
    while (i > 0 && s->dist [i] < s->dist [i - 1]) {
        int ti = s->idx [i];
        double td = s->dist [i];
        s->idx [i] = s->idx [i - 1];
        s->dist [i] = s->dist [i - 1];
        s->idx [i - 1] = ti;
        s->dist [i - 1] = td;
        i--;
    }
}
static void ge_knn_search (const ml_kdnode_t *nodes, int node, ge_knn_t *s) {
    int axis;
    double diff, d2;
    if (node < 0) return;
    d2 = ge_dist2 (&s->pts [(size_t) nodes [node].point * (size_t) s->dim], s->q, s->dim);
    ge_knn_insert (s, nodes [node].point, d2);
    axis = nodes [node].axis;
    diff = s->q [axis] - s->pts [(size_t) nodes [node].point * (size_t) s->dim + (size_t) axis];
    if (diff < 0.0) {
        ge_knn_search (nodes, nodes [node].left, s);
        if (s->count < s->k || diff * diff < s->dist [s->count - 1])
            ge_knn_search (nodes, nodes [node].right, s);
    } else {
        ge_knn_search (nodes, nodes [node].right, s);
        if (s->count < s->k || diff * diff < s->dist [s->count - 1])
            ge_knn_search (nodes, nodes [node].left, s);
    }
}
ML_API int ml_kdtree_knn (const double *pts, const ml_kdnode_t *nodes, int root,
                          int dim, const double *q, int k, int *out_idx,
                          double *out_dist) {
    ge_knn_t s;
    int i;
    if (ML_UNLIKELY (!pts || !nodes || !q || !out_idx || !out_dist || k <= 0 || dim <= 0 || dim > 16))
        return 0;
    if (k > 64) k = 64;
    s.pts = pts;
    s.dim = dim;
    s.k = k;
    s.count = 0;
    for (i = 0; i < dim; i++) s.q [i] = q [i];
    ge_knn_search (nodes, root, &s);
    for (i = 0; i < s.count; i++) {
        out_idx [i] = s.idx [i];
        out_dist [i] = s.dist [i];
    }
    return s.count;
}
static void ge_rad_search (const ml_kdnode_t *nodes, int node, const double *pts,
                           int dim, const double *q, double r2, int *out,
                           int *nout, int maxout) {
    double d2, diff;
    int axis;
    if (node < 0) return;
    d2 = ge_dist2 (&pts [(size_t) nodes [node].point * (size_t) dim], q, dim);
    if (d2 <= r2 && *nout < maxout) out [(*nout)++] = nodes [node].point;
    axis = nodes [node].axis;
    diff = q [axis] - pts [(size_t) nodes [node].point * (size_t) dim + (size_t) axis];
    if (diff < 0.0) {
        ge_rad_search (nodes, nodes [node].left, pts, dim, q, r2, out, nout, maxout);
        if (diff * diff <= r2) ge_rad_search (nodes, nodes [node].right, pts, dim, q, r2, out, nout, maxout);
    } else {
        ge_rad_search (nodes, nodes [node].right, pts, dim, q, r2, out, nout, maxout);
        if (diff * diff <= r2) ge_rad_search (nodes, nodes [node].left, pts, dim, q, r2, out, nout, maxout);
    }
}
ML_API int ml_kdtree_radius (const double *pts, const ml_kdnode_t *nodes,
                             int root, int dim, const double *q, double radius,
                             int *out_idx, int max_out) {
    int nout = 0;
    if (ML_UNLIKELY (!pts || !nodes || !q || !out_idx || dim <= 0 || radius < 0.0))
        return 0;
    ge_rad_search (nodes, root, pts, dim, q, radius * radius, out_idx, &nout, max_out);
    return nout;
}
/* ---------------------------------------------------------- curves ------- */

ML_API void ml_bezier_eval (const double *ctrl, int n, int dim, double t,
                            double *out) {
    double tmp [32 * 8];
    int i, j, d;
    if (ML_UNLIKELY (!ctrl || !out || n <= 0 || n > 32 || dim <= 0 || dim > 8)) return;
    for (i = 0; i < n; i++)
        for (d = 0; d < dim; d++)
            tmp [(size_t) i * (size_t) dim + (size_t) d] = ctrl [(size_t) i * (size_t) dim + (size_t) d];
    for (j = 1; j < n; j++)
        for (i = 0; i < n - j; i++)
            for (d = 0; d < dim; d++)
                tmp [(size_t) i * (size_t) dim + (size_t) d] =
                    (1.0 - t) * tmp [(size_t) i * (size_t) dim + (size_t) d] + t * tmp [(size_t) (i + 1) * (size_t) dim + (size_t) d];
    for (d = 0; d < dim; d++) out [d] = tmp [(size_t) 0 * (size_t) dim + (size_t) d];
}
ML_API void ml_bezier_deriv (const double *ctrl, int n, int dim, double t,
                             double *out) {
    double dc [32 * 8];
    int i, d;
    if (ML_UNLIKELY (!ctrl || !out || n <= 1 || n > 32 || dim <= 0 || dim > 8)) return;
    for (i = 0; i < n - 1; i++)
        for (d = 0; d < dim; d++)
            dc [(size_t) i * (size_t) dim + (size_t) d] =
                (double) (n - 1) * (ctrl [(size_t) (i + 1) * (size_t) dim + (size_t) d] - ctrl [(size_t) i * (size_t) dim + (size_t) d]);
    ml_bezier_eval (dc, n - 1, dim, t, out);
}
ML_API void ml_bspline_eval (const double *ctrl, int nctrl, int degree,
                             const double *knots, double u, int dim,
                             double *out) {
    double d [16 * 8];
    int span, i, j, r, dd;
    int nk = nctrl + degree + 1;
    if (ML_UNLIKELY (!ctrl || !knots || !out || nctrl <= 0 || degree <= 0 || degree > 15 || dim <= 0 || dim > 8))
        return;
    if (u <= knots [degree]) span = degree;
    else if (u >= knots [nctrl]) span = nctrl - 1;
    else {
        span = degree;
        for (i = degree; i < nctrl; i++)
            if (u < knots [i + 1]) {
                span = i;
                break;
            }
    }
    for (j = 0; j <= degree; j++)
        for (dd = 0; dd < dim; dd++)
            d [(size_t) j * (size_t) dim + (size_t) dd] =
                ctrl [(size_t) (span - degree + j) * (size_t) dim + (size_t) dd];
    for (r = 1; r <= degree; r++)
        for (j = degree; j >= r; j--) {
            int kj = span - degree + j;
            double den = knots [kj + degree - r + 1] - knots [kj];
            double a = (den != 0.0) ? (u - knots [kj]) / den : 0.0;
            for (dd = 0; dd < dim; dd++)
                d [(size_t) j * (size_t) dim + (size_t) dd] =
                    (1.0 - a) * d [(size_t) (j - 1) * (size_t) dim + (size_t) dd] + a * d [(size_t) j * (size_t) dim + (size_t) dd];
        }
    for (dd = 0; dd < dim; dd++) out [dd] = d [(size_t) degree * (size_t) dim + (size_t) dd];
    (void) nk;
}
ML_API void ml_catmull_rom (const double *p0, const double *p1,
                            const double *p2, const double *p3, double t,
                            int dim, double *out) {
    double t2 = t * t, t3 = t2 * t;
    int d;
    if (ML_UNLIKELY (!p0 || !p1 || !p2 || !p3 || !out || dim <= 0)) return;
    for (d = 0; d < dim; d++) {
        double a = -0.5 * p0 [d] + 1.5 * p1 [d] - 1.5 * p2 [d] + 0.5 * p3 [d];
        double b = p0 [d] - 2.5 * p1 [d] + 2.0 * p2 [d] - 0.5 * p3 [d];
        double c = -0.5 * p0 [d] + 0.5 * p2 [d];
        out [d] = a * t3 + b * t2 + c * t + p1 [d];
    }
}
/* ---------------------------------------------------------- conics ------- */

ML_API ml_conic_type_t ml_conic_classify (double A, double B, double C,
                                          double D, double E, double F) {
    double disc = B * B - 4.0 * A * C;
    double det;
    if (ml_fabs (A) < 1e-12 && ml_fabs (B) < 1e-12 && ml_fabs (C) < 1e-12) {
        if (ml_fabs (D) < 1e-12 && ml_fabs (E) < 1e-12) return ML_CONIC_EMPTY;
        return ML_CONIC_LINE;
    }
    det = A * (C * F - 0.25 * E * E) - 0.5 * B * (0.5 * B * F - 0.25 * E * D) + 0.5 * D * (0.25 * B * E - 0.5 * C * D);
    if (ml_fabs (disc) < 1e-9) {
        if (ml_fabs (det) < 1e-9) return ML_CONIC_LINE;
        return ML_CONIC_PARABOLA;
    }
    if (disc < 0.0) {
        if (ml_fabs (det) < 1e-9) {
            /* degenerate ellipse: a point or empty */
            double ctrx = (B * E - 2.0 * C * D) / disc;
            double ctry = (B * D - 2.0 * A * E) / disc;
            double val = A * ctrx * ctrx + B * ctrx * ctry + C * ctry * ctry + D * ctrx + E * ctry + F;
            return (val < 0.0) ? ML_CONIC_POINT : ML_CONIC_EMPTY;
        }
        if (ml_fabs (A - C) < 1e-9 && ml_fabs (B) < 1e-9) return ML_CONIC_CIRCLE;
        return ML_CONIC_ELLIPSE;
    }
    if (ml_fabs (det) < 1e-9) return ML_CONIC_LINE;
    return ML_CONIC_HYPERBOLA;
}
/* ------------------------------------------------------- polygons -------- */

ML_API int ml_polygon_clip (const double *sx, const double *sy, int n,
                            const double *cx, const double *cy, int m,
                            double *ox, double *oy, int max_out) {
    double px [512], py [512], qx [512], qy [512];
    int np, i, e;
    if (ML_UNLIKELY (!sx || !sy || !cx || !cy || !ox || !oy || n <= 0 || m < 3 || n > 512))
        return 0;
    np = n;
    for (i = 0; i < n; i++) {
        px [i] = sx [i];
        py [i] = sy [i];
    }
    for (e = 0; e < m; e++) {
        double ax = cx [e], ay = cy [e];
        double bx = cx [(e + 1) % m], by = cy [(e + 1) % m];
        int out = 0;
        for (i = 0; i < np; i++) {
            int k = (i + 1) % np;
            double s = (bx - ax) * (py [i] - ay) - (by - ay) * (px [i] - ax);
            double sk = (bx - ax) * (py [k] - ay) - (by - ay) * (px [k] - ax);
            if (s >= 0.0) {
                qx [out] = px [i];
                qy [out] = py [i];
                out++;
            }
            if ((s > 0.0 && sk < 0.0) || (s < 0.0 && sk > 0.0)) {
                double tt = s / (s - sk);
                qx [out] = px [i] + tt * (px [k] - px [i]);
                qy [out] = py [i] + tt * (py [k] - py [i]);
                out++;
            }
        }
        np = out;
        for (i = 0; i < np; i++) {
            px [i] = qx [i];
            py [i] = qy [i];
        }
        if (np == 0) return 0;
    }
    if (np > max_out) np = max_out;
    for (i = 0; i < np; i++) {
        ox [i] = px [i];
        oy [i] = py [i];
    }
    return np;
}
/* ------------------------------------------------------------ 3-D -------- */

ML_API void ml_seg_seg_closest (const double *a0, const double *a1,
                                const double *b0, const double *b1,
                                double *pa, double *pb) {
    double d1 [3], d2 [3], r [3], a, e, f, s, t, c, b;
    double denom;
    int i;
    for (i = 0; i < 3; i++) {
        d1 [i] = a1 [i] - a0 [i];
        d2 [i] = b1 [i] - b0 [i];
        r [i] = a0 [i] - b0 [i];
    }
    a = d1 [0] * d1 [0] + d1 [1] * d1 [1] + d1 [2] * d1 [2];
    e = d2 [0] * d2 [0] + d2 [1] * d2 [1] + d2 [2] * d2 [2];
    f = d2 [0] * r [0] + d2 [1] * r [1] + d2 [2] * r [2];
    c = d1 [0] * r [0] + d1 [1] * r [1] + d1 [2] * r [2];
    b = d1 [0] * d2 [0] + d1 [1] * d2 [1] + d1 [2] * d2 [2];
    denom = a * e - b * b;
    if (denom > 1e-15) s = (b * f - c * e) / denom;
    else s = 0.0;
    if (s < 0.0) s = 0.0;
    else if (s > 1.0) s = 1.0;
    t = (b * s + f) / (e > 1e-15 ? e : 1e-15);
    if (t < 0.0) {
        t = 0.0;
        s = (-c) / (a > 1e-15 ? a : 1e-15);
        if (s < 0.0) s = 0.0;
        else if (s > 1.0) s = 1.0;
    } else if (t > 1.0) {
        t = 1.0;
        s = (b - c) / (a > 1e-15 ? a : 1e-15);
        if (s < 0.0) s = 0.0;
        else if (s > 1.0) s = 1.0;
    }
    for (i = 0; i < 3; i++) {
        pa [i] = a0 [i] + s * d1 [i];
        pb [i] = b0 [i] + t * d2 [i];
    }
}
