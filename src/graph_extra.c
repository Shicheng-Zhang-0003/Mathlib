#include "ml_compiler.h"
#include "ml_graph_extra.h"
/* ============================================================================
 * v12R2 graph algorithm depth pack.
 * ========================================================================== */

static int ge_dsu_find (int *p, int x) {
    while (p [x] != x) {
        p [x] = p [p [x]];
        x = p [x];
    }
    return x;
}
static int ge_components (int n, const int *eu, const int *ev, int m, int skipv,
                          int skipe) {
    int p [ML_GRAPH_MAXV];
    int i, c;
    if (n > ML_GRAPH_MAXV) return -1;
    for (i = 0; i < n; i++) p [i] = i;
    for (i = 0; i < m; i++) {
        int a, b;
        if (i == skipe) continue;
        if (skipv >= 0 && (eu [i] == skipv || ev [i] == skipv)) continue;
        a = ge_dsu_find (p, eu [i]);
        b = ge_dsu_find (p, ev [i]);
        if (a != b) p [a] = b;
    }
    c = 0;
    for (i = 0; i < n; i++)
        if ((skipv < 0 || i != skipv) && ge_dsu_find (p, i) == i) c++;
    return c;
}
/* ------------------------------------------------------------ all-pairs -- */

ML_API ml_status_t ml_floyd_warshall (int n, const double *W, double *dist,
                                      int *next) {
    int i, j, k;
    if (ML_UNLIKELY (!W || !dist || n <= 0)) return ML_ERR_INVALID_ARG;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++) {
            dist [(size_t) i * (size_t) n + (size_t) j] = W [(size_t) i * (size_t) n + (size_t) j];
            if (next) next [(size_t) i * (size_t) n + (size_t) j] =
                          (i != j && W [(size_t) i * (size_t) n + (size_t) j] < ML_GRAPH_INF) ? j : -1;
        }
    for (k = 0; k < n; k++)
        for (i = 0; i < n; i++) {
            if (dist [(size_t) i * (size_t) n + (size_t) k] >= ML_GRAPH_INF) continue;
            for (j = 0; j < n; j++) {
                double d = dist [(size_t) i * (size_t) n + (size_t) k] + dist [(size_t) k * (size_t) n + (size_t) j];
                if (d < dist [(size_t) i * (size_t) n + (size_t) j]) {
                    dist [(size_t) i * (size_t) n + (size_t) j] = d;
                    if (next) next [(size_t) i * (size_t) n + (size_t) j] =
                                  next [(size_t) i * (size_t) n + (size_t) k];
                }
            }
        }
    for (i = 0; i < n; i++)
        if (dist [(size_t) i * (size_t) n + (size_t) i] < 0.0) return ML_ERR_CONVERGENCE;
    return ML_SUCCESS;
}
ML_API ml_status_t ml_transitive_closure (int n, const double *adj,
                                          unsigned char *reach) {
    int i, j, k;
    if (ML_UNLIKELY (!adj || !reach || n <= 0)) return ML_ERR_INVALID_ARG;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            reach [(size_t) i * (size_t) n + (size_t) j] =
                (unsigned char) ((i == j) || adj [(size_t) i * (size_t) n + (size_t) j] != 0.0);
    for (k = 0; k < n; k++)
        for (i = 0; i < n; i++)
            if (reach [(size_t) i * (size_t) n + (size_t) k])
                for (j = 0; j < n; j++)
                    if (reach [(size_t) k * (size_t) n + (size_t) j])
                        reach [(size_t) i * (size_t) n + (size_t) j] = 1;
    return ML_SUCCESS;
}
/* --------------------------------------------------- single source ------- */

ML_API ml_status_t ml_bellman_ford (int n, const int *eu, const int *ev,
                                    const double *w, int m, int src,
                                    double *dist) {
    int i, it;
    if (ML_UNLIKELY (!eu || !ev || !w || !dist || n <= 0 || src < 0 || src >= n))
        return ML_ERR_INVALID_ARG;
    for (i = 0; i < n; i++) dist [i] = (i == src) ? 0.0 : ML_GRAPH_INF;
    for (it = 0; it < n - 1; it++) {
        int changed = 0;
        for (i = 0; i < m; i++) {
            if (dist [eu [i]] + w [i] < dist [ev [i]]) {
                dist [ev [i]] = dist [eu [i]] + w [i];
                changed = 1;
            }
        }
        if (!changed) break;
    }
    for (i = 0; i < m; i++)
        if (dist [eu [i]] + w [i] < dist [ev [i]]) return ML_ERR_CONVERGENCE;
    return ML_SUCCESS;
}
ML_API ml_status_t ml_astar (int n, const int *eu, const int *ev,
                             const double *w, int m, int src, int dst,
                             ml_heur_t h, void *ctx, double *dist) {
    unsigned char closed [ML_GRAPH_MAXV];
    double f [ML_GRAPH_MAXV];
    int i, iter;
    if (ML_UNLIKELY (!eu || !ev || !w || !dist || n <= 0 || n > ML_GRAPH_MAXV))
        return ML_ERR_INVALID_ARG;
    if (src < 0 || src >= n || dst < 0 || dst >= n) return ML_ERR_INVALID_ARG;
    for (i = 0; i < n; i++) {
        dist [i] = ML_GRAPH_INF;
        closed [i] = 0;
        f [i] = ML_GRAPH_INF;
    }
    dist [src] = 0.0;
    f [src] = h ? h (src, dst, ctx) : 0.0;
    for (iter = 0; iter < n; iter++) {
        int u = -1;
        double best = ML_GRAPH_INF;
        for (i = 0; i < n; i++)
            if (!closed [i] && f [i] < best) {
                best = f [i];
                u = i;
            }
        if (u < 0) break;
        closed [u] = 1;
        if (u == dst) break;
        for (i = 0; i < m; i++) {
            int v = -1;
            if (eu [i] == u) v = ev [i];
            else if (ev [i] == u) v = eu [i];
            if (v < 0 || closed [v]) continue;
            if (dist [u] + w [i] < dist [v]) {
                dist [v] = dist [u] + w [i];
                f [v] = dist [v] + (h ? h (v, dst, ctx) : 0.0);
            }
        }
    }
    return ML_SUCCESS;
}
/* ----------------------------------------------------------- max-flow ---- */

ML_API double ml_max_flow (int n, const double *cap, int s, int t,
                           double *flow) {
    double res [64][64];
    int i, j;
    if (ML_UNLIKELY (!cap || n <= 0 || n > 64 || s < 0 || t < 0 || s >= n || t >= n || s == t))
        return 0.0;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++) res [i][j] = cap [(size_t) i * (size_t) n + (size_t) j];
    for (;;) {
        int par [64];
        unsigned char seen [64];
        int q [64], qh = 0, qt = 0, u;
        for (i = 0; i < n; i++) {
            par [i] = -1;
            seen [i] = 0;
        }
        seen [s] = 1;
        q [qt++] = s;
        while (qh < qt && !seen [t]) {
            u = q [qh++];
            for (i = 0; i < n; i++)
                if (!seen [i] && res [u][i] > 1e-12) {
                    seen [i] = 1;
                    par [i] = u;
                    q [qt++] = i;
                    if (i == t) break;
                }
        }
        if (!seen [t]) break;
        {
            double aug = ML_GRAPH_INF;
            for (u = t; u != s; u = par [u]) {
                double r = res [par [u]][u];
                if (r < aug) aug = r;
            }
            for (u = t; u != s; u = par [u]) {
                res [par [u]][u] -= aug;
                res [u][par [u]] += aug;
            }
        }
    }
    if (flow) {
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++) {
                double fv = cap [(size_t) i * (size_t) n + (size_t) j] - res [i][j];
                flow [(size_t) i * (size_t) n + (size_t) j] = fv > 0.0 ? fv : 0.0;
            }
    }
    {
        double mf = 0.0;
        for (i = 0; i < n; i++)
            mf += cap [(size_t) i * (size_t) n + (size_t) t] - res [i][t];
        return mf;
    }
}
/* --------------------------------------------------- bipartite matching -- */

static int ge_try_kuhn (int u, int nR, const int *adj, int *matchR,
                        unsigned char *vis) {
    int v;
    for (v = 0; v < nR; v++) {
        if (!adj [(size_t) u * (size_t) nR + (size_t) v] || vis [v]) continue;
        vis [v] = 1;
        if (matchR [v] < 0 || ge_try_kuhn (matchR [v], nR, adj, matchR, vis)) {
            matchR [v] = u;
            return 1;
        }
    }
    return 0;
}
ML_API int ml_bipartite_matching (int nL, int nR, const int *adj, int *matchL,
                                  int *matchR) {
    int u, res = 0;
    if (ML_UNLIKELY (!adj || nL <= 0 || nR <= 0)) return 0;
    for (u = 0; u < nR; u++) matchR [u] = -1;
    for (u = 0; u < nL; u++) {
        unsigned char vis [ML_GRAPH_MAXV];
        int i;
        for (i = 0; i < nR; i++) vis [i] = 0;
        if (ge_try_kuhn (u, nR, adj, matchR, vis)) res++;
    }
    if (matchL)
        for (u = 0; u < nL; u++) matchL [u] = -1;
    if (matchL)
        for (u = 0; u < nR; u++)
            if (matchR [u] >= 0) matchL [matchR [u]] = u;
    return res;
}
/* -------------------------------------------- connectivity structure ----- */

static void ge_tarjan (int n, const int *eu, const int *ev, int m, int u,
                       int *idx, int *low, int *onstk, int *stk, int *sp,
                       int *comp, int *cid, int *timer) {
    int k;
    idx [u] = low [u] = (*timer)++;
    stk [(*sp)++] = u;
    onstk [u] = 1;
    for (k = 0; k < m; k++) {
        int v;
        if (eu [k] != u) continue;
        v = ev [k];
        if (idx [v] < 0) {
            ge_tarjan (n, eu, ev, m, v, idx, low, onstk, stk, sp, comp, cid, timer);
            if (low [v] < low [u]) low [u] = low [v];
        } else if (onstk [v]) {
            if (idx [v] < low [u]) low [u] = idx [v];
        }
    }
    if (low [u] == idx [u]) {
        int w;
        do {
            w = stk [--(*sp)];
            onstk [w] = 0;
            comp [w] = *cid;
        } while (w != u);
        (*cid)++;
    }
}
ML_API int ml_scc (int n, const int *eu, const int *ev, int m, int *comp) {
    int idx [ML_GRAPH_MAXV], low [ML_GRAPH_MAXV], onstk [ML_GRAPH_MAXV];
    int stk [ML_GRAPH_MAXV], sp = 0, cid = 0, timer = 0, u;
    if (ML_UNLIKELY (!eu || !ev || !comp || n <= 0 || n > ML_GRAPH_MAXV))
        return 0;
    for (u = 0; u < n; u++) {
        idx [u] = -1;
        onstk [u] = 0;
        comp [u] = -1;
    }
    for (u = 0; u < n; u++)
        if (idx [u] < 0)
            ge_tarjan (n, eu, ev, m, u, idx, low, onstk, stk, &sp, comp, &cid, &timer);
    return cid;
}
ML_API ml_status_t ml_articulation_points (int n, const int *eu, const int *ev,
                                           int m, unsigned char *is_ap,
                                           unsigned char *is_bridge) {
    int base, v, e;
    if (ML_UNLIKELY (!eu || !ev || n <= 0 || n > ML_GRAPH_MAXV))
        return ML_ERR_INVALID_ARG;
    base = ge_components (n, eu, ev, m, -1, -1);
    if (is_ap)
        for (v = 0; v < n; v++) {
            int after = ge_components (n, eu, ev, m, v, -1);
            is_ap [v] = (unsigned char) (after > base);
        }
    if (is_bridge)
        for (e = 0; e < m; e++) {
            int after = ge_components (n, eu, ev, m, -1, e);
            is_bridge [e] = (unsigned char) (after > base);
        }
    return ML_SUCCESS;
}
/* --------------------------------------------------------- Euler path ---- */

ML_API int ml_eulerian_path (int n, const int *eu, const int *ev, int m,
                             int *path) {
    unsigned char used [4096];
    int deg [ML_GRAPH_MAXV], stk [ML_GRAPH_MAXV + 4096], sp = 0;
    int odd = 0, start = -1, i, out = 0, nz = 0;
    if (ML_UNLIKELY (!eu || !ev || !path || n <= 0 || n > ML_GRAPH_MAXV || m < 0 || m > 4096))
        return 0;
    for (i = 0; i < n; i++) deg [i] = 0;
    for (i = 0; i < m; i++) {
        deg [eu [i]]++;
        deg [ev [i]]++;
        used [i] = 0;
    }
    for (i = 0; i < n; i++) {
        if (deg [i] & 1) {
            odd++;
            start = i;
        }
        if (deg [i] > 0) nz++;
    }
    if (odd != 0 && odd != 2) return 0;
    if (start < 0) {
        for (i = 0; i < n; i++)
            if (deg [i] > 0) {
                start = i;
                break;
            }
    }
    if (start < 0) return (m == 0) ? 0 : 0;
    /* Reachability of all nonzero-degree vertices from start. */
    {
        unsigned char seen [ML_GRAPH_MAXV];
        int st [ML_GRAPH_MAXV], ss = 0, reached = 0;
        for (i = 0; i < n; i++) seen [i] = 0;
        seen [start] = 1;
        st [ss++] = start;
        while (ss > 0) {
            int v = st [--ss], k;
            for (k = 0; k < m; k++) {
                int w = -1;
                if (eu [k] == v) w = ev [k];
                else if (ev [k] == v) w = eu [k];
                if (w >= 0 && !seen [w]) {
                    seen [w] = 1;
                    st [ss++] = w;
                }
            }
        }
        for (i = 0; i < n; i++)
            if (deg [i] > 0 && seen [i]) reached++;
        if (reached != nz) return 0;
    }
    stk [sp++] = start;
    while (sp > 0) {
        int v = stk [sp - 1], found = -1, other = -1;
        for (i = 0; i < m; i++) {
            if (used [i]) continue;
            if (eu [i] == v) {
                found = i;
                other = ev [i];
                break;
            }
            if (ev [i] == v) {
                found = i;
                other = eu [i];
                break;
            }
        }
        if (found >= 0) {
            used [found] = 1;
            stk [sp++] = other;
        } else {
            path [out++] = v;
            sp--;
        }
    }
    if (out != m + 1) return 0;
    for (i = 0; i < out / 2; i++) {
        int t = path [i];
        path [i] = path [out - 1 - i];
        path [out - 1 - i] = t;
    }
    return out;
}
/* ------------------------------------------------------------ PageRank --- */

ML_API ml_status_t ml_pagerank (int n, const int *eu, const int *ev, int m,
                                double damping, int iters, double *pr) {
    double out [ML_GRAPH_MAXV], nxt [ML_GRAPH_MAXV];
    int i, it;
    if (ML_UNLIKELY (!eu || !ev || !pr || n <= 0 || n > ML_GRAPH_MAXV))
        return ML_ERR_INVALID_ARG;
    if (damping <= 0.0 || damping >= 1.0) damping = 0.85;
    if (iters <= 0) iters = 100;
    for (i = 0; i < n; i++) {
        out [i] = 0.0;
        pr [i] = 1.0 / (double) n;
    }
    for (i = 0; i < m; i++) out [eu [i]] += 1.0;
    for (it = 0; it < iters; it++) {
        double dangling = 0.0, base;
        for (i = 0; i < n; i++) {
            nxt [i] = 0.0;
            if (out [i] == 0.0) dangling += pr [i];
        }
        base = (1.0 - damping) / (double) n + damping * dangling / (double) n;
        for (i = 0; i < n; i++) nxt [i] = base;
        for (i = 0; i < m; i++)
            nxt [ev [i]] += damping * pr [eu [i]] / out [eu [i]];
        for (i = 0; i < n; i++) pr [i] = nxt [i];
    }
    return ML_SUCCESS;
}
