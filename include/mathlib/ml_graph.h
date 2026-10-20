#ifndef MATHLIB_ML_GRAPH_H
#define MATHLIB_ML_GRAPH_H

#include "ml_compiler.h"
#include "ml_core.h"
#include "ml_types.h"
#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * Graph theory: union-find, traversal, shortest paths, spanning trees.
 *
 * Representation is edge lists (eu/ev endpoint arrays, optional double
 * weights); the library builds what it needs per call. Conventions:
 * - Vertices are 0..n-1; endpoints outside range -> ML_ERR_INVALID_ARG.
 * - Caps: n in 1..4096, m in 0..262144 (heap scratch with SIZE_MAX guards;
 *   Dijkstra is O(n^2) simple-scan, Kruskal sorts via qsort).
 * - Self-loops are ignored by traversal/MST (documented); parallel edges
 *   are kept (minimum wins in Dijkstra/Kruskal naturally).
 * - Dijkstra/Kruskal need finite weights >= 0 (negative -> INVALID_ARG;
 *   Bellman-Ford is deferred, documented).
 * - int verdicts fail closed (0 / -1); dist arrays use -1.0 for
 *   unreachable (all weights >= 0, so -1.0 is unambiguous).
 * ========================================================================== */

/* Disjoint set union (caller-owned arrays, no shared state) */
ML_API void ml_dsu_make (int *parent, int *rank, int n);
ML_API int ml_dsu_find (int *parent, int n, int x);
ML_API void ml_dsu_union (int *parent, int *rank, int n, int a, int b);

/* Traversal + structure (edge lists, undirected unless noted) */
ML_API ml_status_t ml_bfs_dist (int n, const int *eu, const int *ev, int m,
                                int src, double *dist);
ML_API ml_status_t ml_dfs_order (int n, const int *eu, const int *ev, int m,
                                 int src, int *order, int *count);
ML_API int ml_is_bipartite (int n, const int *eu, const int *ev, int m);
ML_API int ml_components (int n, const int *eu, const int *ev, int m);
/* Kahn topological sort on DIRECTED edges (u -> v); returns count, or -1
 * when a cycle exists (out untouched on cycle). Self-loop = cycle. */
ML_API int ml_topo_sort (int n, const int *eu, const int *ev, int m, int *out);

/* Weighted (w[i] >= 0 finite per edge i). Kruskal succeeds only on
 * connected graphs (SINGULAR when disconnected). */
ML_API ml_status_t ml_dijkstra (int n, const int *eu, const int *ev,
                                const double *w, int m, int src, double *dist);
ML_API ml_status_t ml_kruskal_mst (int n, const int *eu, const int *ev,
                                   const double *w, int m, double *total);

#ifdef __cplusplus
}
#endif
#endif /* MATHLIB_ML_GRAPH_H */
