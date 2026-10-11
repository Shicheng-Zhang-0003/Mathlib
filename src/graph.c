#include "ml_compiler.h"
#include "ml_graph.h"
#include "ml_core.h"
#include <stddef.h>
#include <stdlib.h>

/* Graph theory on edge lists. Caller owns all arrays; scratch is
 * heap-per-call with SIZE_MAX guards (no shared mutable state). */

#define MLG_NMAX 4096
#define MLG_MMAX 262144

static int mlg_valid_edge_args(int n, const int *eu, const int *ev, int m) {
    if (!eu || !ev || n < 1 || n > MLG_NMAX || m < 0 || m > MLG_MMAX) return 0;
    return 1;
}

static int mlg_check_endpoints(int n, const int *eu, const int *ev, int m) {
    int i;
    for (i = 0; i < m; i++) {
        if (eu[i] < 0 || eu[i] >= n || ev[i] < 0 || ev[i] >= n) return 0;
    }
    return 1;
}

/* CSR builder: head/to/nxt over 2*m directed arcs (self-loops kept once).
 * Returns 0 on OOM (caller maps to status). */
static int mlg_build_csr(int n, const int *eu, const int *ev, int m,
                         int **head, int **to, int **nxt, int *narc) {
    int i;
    int e = 0;
    size_t nn = (size_t)n;
    if (nn > (size_t)-1 / sizeof(int)) return 0;
    if ((size_t)m > (size_t)-1 / (2 * sizeof(int))) return 0;
    *head = (int *)malloc(nn * sizeof(int));
    *to = (int *)malloc((size_t)(2 * m) * sizeof(int));
    *nxt = (int *)malloc((size_t)(2 * m) * sizeof(int));
    if (!*head || ((m != 0) && (!*to || !*nxt))) {
        free(*head); free(*to); free(*nxt);
        *head = NULL; *to = NULL; *nxt = NULL;
        return 0;
    }
    for (i = 0; i < n; i++) (*head)[i] = -1;
    for (i = 0; i < m; i++) {
        int u = eu[i], v = ev[i];
        (*to)[e] = v; (*nxt)[e] = (*head)[u]; (*head)[u] = e; e++;
        if (v != u) {
            (*to)[e] = u; (*nxt)[e] = (*head)[v]; (*head)[v] = e; e++;
        }
    }
    *narc = e;
    return 1;
}

ML_API void ml_dsu_make(int *parent, int *rank, int n) {
    int i;
    if (!parent || !rank || n <= 0) return;
    for (i = 0; i < n; i++) { parent[i] = i; rank[i] = 0; }
}

ML_API int ml_dsu_find(int *parent, int n, int x) {
    int root;
    if (!parent || n <= 0 || x < 0 || x >= n) return -1;
    root = x;
    while (parent[root] != root) {
        if (parent[root] < 0 || parent[root] >= n) return -1;
        root = parent[root];
    }
    while (parent[x] != root) {
        int nxt = parent[x];
        parent[x] = root;
        x = nxt;
    }
    return root;
}

ML_API void ml_dsu_union(int *parent, int *rank, int n, int a, int b) {
    int ra;
    int rb;
    if (!parent || !rank || n <= 0) return;
    ra = ml_dsu_find(parent, n, a);
    rb = ml_dsu_find(parent, n, b);
    if (ra < 0 || rb < 0 || ra == rb) return;
    if (rank[ra] < rank[rb]) parent[ra] = rb;
    else if (rank[ra] > rank[rb]) parent[rb] = ra;
    else { parent[rb] = ra; rank[ra]++; }
}

ML_API ml_status_t ml_bfs_dist(int n, const int *eu, const int *ev, int m,
                               int src, double *dist) {
    int *head = NULL, *to = NULL, *nxt = NULL;
    int *queue = NULL;
    int narc = 0;
    int qh = 0, qt = 0;
    int i;
    if (!mlg_valid_edge_args(n, eu, ev, m)) return ML_ERR_INVALID_ARG;
    if (!dist || src < 0 || src >= n) return ML_ERR_INVALID_ARG;
    if (!mlg_check_endpoints(n, eu, ev, m)) return ML_ERR_INVALID_ARG;
    if ((size_t)n > (size_t)-1 / sizeof(int)) return ML_ERR_INVALID_ARG;
    queue = (int *)malloc((size_t)n * sizeof(int));
    if (!queue) return ML_ERR_WORKSPACE;
    if (!mlg_build_csr(n, eu, ev, m, &head, &to, &nxt, &narc)) {
        free(queue);
        return ML_ERR_WORKSPACE;
    }
    for (i = 0; i < n; i++) dist[i] = -1.0;
    dist[src] = 0.0;
    queue[qt++] = src;
    while (qh < qt) {
        int u = queue[qh++];
        int e;
        for (e = head[u]; e >= 0; e = nxt[e]) {
            int v = to[e];
            if (dist[v] < 0.0) {
                dist[v] = dist[u] + 1.0;
                queue[qt++] = v;
            }
        }
    }
    free(queue);
    free(head); free(to); free(nxt);
    return ML_SUCCESS;
}

ML_API ml_status_t ml_dfs_order(int n, const int *eu, const int *ev, int m,
                                int src, int *order, int *count) {
    int *head = NULL, *to = NULL, *nxt = NULL;
    int *stack = NULL;
    int *seen = NULL;
    int narc = 0;
    int top = 0;
    int c = 0;
    int i;
    if (!mlg_valid_edge_args(n, eu, ev, m)) return ML_ERR_INVALID_ARG;
    if (!order || !count || src < 0 || src >= n) return ML_ERR_INVALID_ARG;
    if (!mlg_check_endpoints(n, eu, ev, m)) return ML_ERR_INVALID_ARG;
    if ((size_t)n > (size_t)-1 / sizeof(int)) return ML_ERR_INVALID_ARG;
    stack = (int *)malloc((size_t)n * sizeof(int));
    seen = (int *)malloc((size_t)n * sizeof(int));
    if (!stack || !seen) {
        free(stack); free(seen);
        return ML_ERR_WORKSPACE;
    }
    if (!mlg_build_csr(n, eu, ev, m, &head, &to, &nxt, &narc)) {
        free(stack); free(seen);
        return ML_ERR_WORKSPACE;
    }
    for (i = 0; i < n; i++) seen[i] = 0;
    stack[top++] = src;
    seen[src] = 1;
    while (top > 0) {
        int u = stack[--top];
        int e;
        order[c++] = u;
        for (e = head[u]; e >= 0; e = nxt[e]) {
            int v = to[e];
            if (!seen[v]) {
                seen[v] = 1;
                stack[top++] = v;
            }
        }
    }
    *count = c;
    free(stack); free(seen);
    free(head); free(to); free(nxt);
    return ML_SUCCESS;
}

ML_API int ml_is_bipartite(int n, const int *eu, const int *ev, int m) {
    int *head = NULL, *to = NULL, *nxt = NULL;
    int *color = NULL;
    int *queue = NULL;
    int narc = 0;
    int i;
    int ok = 1;
    if (!mlg_valid_edge_args(n, eu, ev, m)) return 0;
    if (!mlg_check_endpoints(n, eu, ev, m)) return 0;
    if ((size_t)n > (size_t)-1 / sizeof(int)) return 0;
    color = (int *)malloc((size_t)n * sizeof(int));
    queue = (int *)malloc((size_t)n * sizeof(int));
    if (!color || !queue) {
        free(color); free(queue);
        return 0;
    }
    if (!mlg_build_csr(n, eu, ev, m, &head, &to, &nxt, &narc)) {
        free(color); free(queue);
        return 0;
    }
    for (i = 0; i < n; i++) color[i] = -1;
    for (i = 0; i < n && ok; i++) {
        int qh = 0, qt = 0;
        if (color[i] != -1) continue;
        color[i] = 0;
        queue[qt++] = i;
        while (qh < qt && ok) {
            int u = queue[qh++];
            int e;
            for (e = head[u]; e >= 0; e = nxt[e]) {
                int v = to[e];
                if (v == u) { ok = 0; break; }
                if (color[v] == -1) {
                    color[v] = color[u] ^ 1;
                    queue[qt++] = v;
                } else if (color[v] == color[u]) {
                    ok = 0;
                    break;
                }
            }
        }
    }
    free(color); free(queue);
    free(head); free(to); free(nxt);
    return ok;
}

ML_API int ml_components(int n, const int *eu, const int *ev, int m) {
    int *parent = NULL;
    int *rank = NULL;
    int cnt = 0;
    int i;
    if (!mlg_valid_edge_args(n, eu, ev, m)) return -1;
    if (!mlg_check_endpoints(n, eu, ev, m)) return -1;
    if ((size_t)n > (size_t)-1 / sizeof(int)) return -1;
    parent = (int *)malloc((size_t)n * sizeof(int));
    rank = (int *)malloc((size_t)n * sizeof(int));
    if (!parent || !rank) {
        free(parent); free(rank);
        return -1;
    }
    ml_dsu_make(parent, rank, n);
    for (i = 0; i < m; i++) {
        if (eu[i] != ev[i]) ml_dsu_union(parent, rank, n, eu[i], ev[i]);
    }
    for (i = 0; i < n; i++) {
        if (ml_dsu_find(parent, n, i) == i) cnt++;
    }
    free(parent); free(rank);
    return cnt;
}

ML_API int ml_topo_sort(int n, const int *eu, const int *ev, int m, int *out) {
    int *indeg = NULL;
    int *queue = NULL;
    int qh = 0, qt = 0;
    int c = 0;
    int i;
    if (!mlg_valid_edge_args(n, eu, ev, m)) return -1;
    if (!out) return -1;
    if (!mlg_check_endpoints(n, eu, ev, m)) return -1;
    if ((size_t)n > (size_t)-1 / sizeof(int)) return -1;
    indeg = (int *)malloc((size_t)n * sizeof(int));
    queue = (int *)malloc((size_t)n * sizeof(int));
    if (!indeg || !queue) {
        free(indeg); free(queue);
        return -1;
    }
    for (i = 0; i < n; i++) indeg[i] = 0;
    for (i = 0; i < m; i++) indeg[ev[i]]++;
    for (i = 0; i < n; i++) {
        if (indeg[i] == 0) queue[qt++] = i;
    }
    while (qh < qt) {
        int u = queue[qh++];
        out[c++] = u;
        for (i = 0; i < m; i++) {
            if (eu[i] == u) {
                if (--indeg[ev[i]] == 0) queue[qt++] = ev[i];
            }
        }
    }
    free(indeg); free(queue);
    return (c == n) ? c : -1;
}

static int mlg_check_weights(const double *w, int m) {
    int i;
    if (!w) return 0;
    for (i = 0; i < m; i++) {
        if (ml_isnan(w[i]) || !ml_isfinite(w[i]) || w[i] < 0.0) return 0;
    }
    return 1;
}

ML_API ml_status_t ml_dijkstra(int n, const int *eu, const int *ev,
                               const double *w, int m, int src, double *dist) {
    int *head = NULL, *to = NULL, *nxt = NULL;
    double *ew = NULL;
    int *done = NULL;
    int narc = 0;
    int i;
    if (!mlg_valid_edge_args(n, eu, ev, m)) return ML_ERR_INVALID_ARG;
    if (!dist || src < 0 || src >= n) return ML_ERR_INVALID_ARG;
    if (!mlg_check_endpoints(n, eu, ev, m)) return ML_ERR_INVALID_ARG;
    if (!mlg_check_weights(w, m)) return ML_ERR_INVALID_ARG;
    if ((size_t)n > (size_t)-1 / sizeof(int)) return ML_ERR_INVALID_ARG;
    if ((size_t)(2 * m) > (size_t)-1 / sizeof(double) && m > 0) return ML_ERR_INVALID_ARG;
    done = (int *)malloc((size_t)n * sizeof(int));
    ew = (m == 0) ? NULL : (double *)malloc((size_t)(2 * m) * sizeof(double));
    if (!done || (m != 0 && !ew)) {
        free(done); free(ew);
        return ML_ERR_WORKSPACE;
    }
    if (!mlg_build_csr(n, eu, ev, m, &head, &to, &nxt, &narc)) {
        free(done); free(ew);
        return ML_ERR_WORKSPACE;
    }
    /* Mirror edge weights onto arcs: arc order in mlg_build_csr is
     * deterministic (edge i then its reverse), so replay it. */
    {
        int e = 0;
        for (i = 0; i < m; i++) {
            ew[e++] = w[i];
            if (ev[i] != eu[i]) ew[e++] = w[i];
        }
    }
    for (i = 0; i < n; i++) { dist[i] = -1.0; done[i] = 0; }
    dist[src] = 0.0;
    for (i = 0; i < n; i++) {
        int u = -1;
        int j;
        int e;
        double best = -1.0;
        for (j = 0; j < n; j++) {
            if (!done[j] && dist[j] >= 0.0 && (u < 0 || dist[j] < best)) {
                u = j;
                best = dist[j];
            }
        }
        if (u < 0) break;
        done[u] = 1;
        for (e = head[u]; e >= 0; e = nxt[e]) {
            int v = to[e];
            double nd = dist[u] + ew[e];
            if (!ml_isfinite(nd)) continue;
            if (dist[v] < 0.0 || nd < dist[v]) dist[v] = nd;
        }
    }
    free(done); free(ew);
    free(head); free(to); free(nxt);
    return ML_SUCCESS;
}

typedef struct {
    int u;
    int v;
    double w;
} mlg_edge_t;

static int mlg_edge_cmp(const void *pa, const void *pb) {
    const mlg_edge_t *a = (const mlg_edge_t *)pa;
    const mlg_edge_t *b = (const mlg_edge_t *)pb;
    if (a->w < b->w) return -1;
    if (a->w > b->w) return 1;
    return 0;
}

ML_API ml_status_t ml_kruskal_mst(int n, const int *eu, const int *ev,
                                  const double *w, int m, double *total) {
    mlg_edge_t *edges = NULL;
    int *parent = NULL;
    int *rank = NULL;
    long double sum = 0.0L;
    int used = 0;
    int i;
    if (!mlg_valid_edge_args(n, eu, ev, m)) return ML_ERR_INVALID_ARG;
    if (!total) return ML_ERR_INVALID_ARG;
    if (!mlg_check_endpoints(n, eu, ev, m)) return ML_ERR_INVALID_ARG;
    if (!mlg_check_weights(w, m)) return ML_ERR_INVALID_ARG;
    if (m > 0 && (size_t)m > (size_t)-1 / sizeof(mlg_edge_t)) return ML_ERR_INVALID_ARG;
    if ((size_t)n > (size_t)-1 / sizeof(int)) return ML_ERR_INVALID_ARG;
    edges = (m == 0) ? NULL : (mlg_edge_t *)malloc((size_t)m * sizeof(mlg_edge_t));
    parent = (int *)malloc((size_t)n * sizeof(int));
    rank = (int *)malloc((size_t)n * sizeof(int));
    if ((m != 0 && !edges) || !parent || !rank) {
        free(edges); free(parent); free(rank);
        return ML_ERR_WORKSPACE;
    }
    for (i = 0; i < m; i++) {
        if (eu[i] == ev[i]) continue;
        edges[used].u = eu[i];
        edges[used].v = ev[i];
        edges[used].w = w[i];
        used++;
    }
    qsort(edges, (size_t)used, sizeof(mlg_edge_t), mlg_edge_cmp);
    ml_dsu_make(parent, rank, n);
    {
        int taken = 0;
        for (i = 0; i < used; i++) {
            int ru = ml_dsu_find(parent, n, edges[i].u);
            int rv = ml_dsu_find(parent, n, edges[i].v);
            if (ru < 0 || rv < 0 || ru == rv) continue;
            ml_dsu_union(parent, rank, n, edges[i].u, edges[i].v);
            sum += (long double)edges[i].w;
            taken++;
            if (taken == n - 1) break;
        }
        free(edges); free(parent); free(rank);
        if (taken != n - 1) return ML_ERR_SINGULAR;
    }
    *total = (double)sum;
    if (!ml_isfinite(*total)) return ML_ERR_OVERFLOW;
    return ML_SUCCESS;
}
