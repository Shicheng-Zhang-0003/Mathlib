#ifndef MATHLIB_ML_GRAPH_EXTRA_H
#define MATHLIB_ML_GRAPH_EXTRA_H

#include "ml_compiler.h"
#include "ml_core.h"
#include "ml_types.h"
#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * v12R2 GRAPH ALGORITHM DEPTH PACK
 *
 * All-pairs shortest paths, negative-weight shortest paths, A*, max-flow,
 * bipartite matching, strongly connected components, articulation points,
 * bridges, Eulerian paths and PageRank.  Graphs are either edge lists
 * (eu/ev arrays, optionally weighted) or dense adjacency matrices.
 *
 * No dynamic allocation: explicit DFS stacks are bounded and large DP tables
 * are caller-supplied.
 * ========================================================================== */

#define ML_GRAPH_INF 1e300
#define ML_GRAPH_MAXV 1024

/* Heuristic for A*: estimated cost from `node` to `goal`. */
typedef double (*ml_heur_t) (int node, int goal, void *ctx);

/* ------------------------------------------------------------ all-pairs -- */
/* W is n x n; W[i*n+j] = weight of edge i->j or ML_GRAPH_INF if none.
 * dist[i*n+j] receives shortest distances; next[i*n+j] the successor on the
 * path (or -1).  Negative cycles produce ML_ERR_CONVERGENCE. */
ML_API ml_status_t ml_floyd_warshall (int n, const double *W, double *dist,
                                      int *next);

/* Boolean reachability from dense adjacency `adj` (n x n, nonzero = edge). */
ML_API ml_status_t ml_transitive_closure (int n, const double *adj,
                                          unsigned char *reach);

/* --------------------------------------------------- single source ------- */
/* Bellman-Ford over edge list.  Returns ML_ERR_CONVERGENCE on a negative
 * cycle reachable from src. */
ML_API ml_status_t ml_bellman_ford (int n, const int *eu, const int *ev,
                                    const double *w, int m, int src,
                                    double *dist);

/* A* over an undirected weighted edge list.  Distances written to dist[n]. */
ML_API ml_status_t ml_astar (int n, const int *eu, const int *ev,
                             const double *w, int m, int src, int dst,
                             ml_heur_t h, void *ctx, double *dist);

/* ----------------------------------------------------------- max-flow ---- */
/* Edmonds-Karp.  cap is the n x n capacity matrix, flow receives the actual
 * flow (may alias NULL to skip).  Returns the maximum s->t flow.  n<=64. */
ML_API double ml_max_flow (int n, const double *cap, int s, int t,
                           double *flow);

/* --------------------------------------------------- bipartite matching -- */
/* adj is nL x nR (nonzero = edge).  Fills matchL (length nL) and matchR
 * (length nR) with -1 for unmatched.  Returns the matching size. */
ML_API int ml_bipartite_matching (int nL, int nR, const int *adj, int *matchL,
                                  int *matchR);

/* -------------------------------------------- connectivity structure ----- */
/* Tarjan strongly connected components over an edge list.  comp[i] gets a
 * component id (0-based, in reverse topological order).  Returns #components. */
ML_API int ml_scc (int n, const int *eu, const int *ev, int m, int *comp);

/* Articulation points and bridges.  is_ap and is_bridge may be NULL. */
ML_API ml_status_t ml_articulation_points (int n, const int *eu, const int *ev,
                                           int m, unsigned char *is_ap,
                                           unsigned char *is_bridge);

/* --------------------------------------------------------- Euler path ---- */
/* Hierholzer.  path must hold m+1 vertex ids.  Returns the number of vertices
 * written (m+1) or 0 when no Eulerian path/circuit exists. */
ML_API int ml_eulerian_path (int n, const int *eu, const int *ev, int m,
                             int *path);

/* ------------------------------------------------------------ PageRank --- */
ML_API ml_status_t ml_pagerank (int n, const int *eu, const int *ev, int m,
                                double damping, int iters, double *pr);

#ifdef __cplusplus
}
#endif
#endif /* MATHLIB_ML_GRAPH_EXTRA_H */
