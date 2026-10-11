/* Graph pins: DSU, BFS/DFS, bipartite, components, topo, Dijkstra, Kruskal.
 * Ground truth: hand-computed small graphs. */
#include "test_harness.h"
#include "ml_graph.h"
#include "ml_core.h"

int main(void) {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Edge Graph");

    /* DSU. */
    {
        int p[4], r[4];
        ml_dsu_make(p, r, 4);
        ml_dsu_union(p, r, 4, 0, 1);
        ASSERT_TRUE(&ctx, ml_dsu_find(p, 4, 0) == ml_dsu_find(p, 4, 1), "dsu union");
        ASSERT_TRUE(&ctx, ml_dsu_find(p, 4, 0) != ml_dsu_find(p, 4, 2), "dsu separate");
        ASSERT_TRUE(&ctx, ml_dsu_find(p, 4, 9) == -1, "dsu out of range");
        ASSERT_TRUE(&ctx, ml_dsu_find(NULL, 4, 0) == -1, "dsu null");
    }

    /* Path 0-1-2-3, weights 1,2,1. */
    {
        int eu[] = {0, 1, 2}, ev[] = {1, 2, 3};
        double w[] = {1.0, 2.0, 1.0};
        double d[4];
        double tot = 0;
        int o[4], c = 0;
        ASSERT_TRUE(&ctx, ml_bfs_dist(4, eu, ev, 3, 0, d) == ML_SUCCESS, "bfs st");
        ASSERT_NEAR(&ctx, d[0], 0.0, 0.0, "bfs d0");
        ASSERT_NEAR(&ctx, d[3], 3.0, 0.0, "bfs d3");
        ASSERT_TRUE(&ctx, ml_dijkstra(4, eu, ev, w, 3, 0, d) == ML_SUCCESS, "dij st");
        ASSERT_NEAR(&ctx, d[2], 3.0, 1e-12, "dij weighted");
        ASSERT_NEAR(&ctx, d[3], 4.0, 1e-12, "dij end");
        ASSERT_TRUE(&ctx, ml_kruskal_mst(4, eu, ev, w, 3, &tot) == ML_SUCCESS, "kr st");
        ASSERT_NEAR(&ctx, tot, 4.0, 1e-12, "kr total");
        ASSERT_TRUE(&ctx, ml_dfs_order(4, eu, ev, 3, 0, o, &c) == ML_SUCCESS && c == 4, "dfs visits all");
        ASSERT_TRUE(&ctx, ml_kruskal_mst(4, eu, ev, w, 2, &tot) == ML_ERR_SINGULAR, "kr disconnected");
    }

    /* Triangle MST picks two cheapest. */
    {
        int eu[] = {0, 1, 1, 0}, ev[] = {1, 2, 2, 2};
        double w[] = {3.0, 1.0, 1.0, 2.0};
        double tot = 0;
        ASSERT_TRUE(&ctx, ml_kruskal_mst(3, eu, ev, w, 4, &tot) == ML_SUCCESS, "kr tri");
        ASSERT_NEAR(&ctx, tot, 3.0, 1e-12, "kr tri total");
    }

    /* Bipartite: even cycle yes, odd cycle no, self-loop no. */
    {
        int e4[] = {0, 1, 2, 3}, f4[] = {1, 2, 3, 0};
        int e3[] = {0, 1, 2}, f3[] = {1, 2, 0};
        int es[] = {0}, fs[] = {0};
        ASSERT_TRUE(&ctx, ml_is_bipartite(4, e4, f4, 4) == 1, "even cycle");
        ASSERT_TRUE(&ctx, ml_is_bipartite(3, e3, f3, 3) == 0, "odd cycle");
        ASSERT_TRUE(&ctx, ml_is_bipartite(1, es, fs, 1) == 0, "self-loop");
        ASSERT_TRUE(&ctx, ml_is_bipartite(3, NULL, NULL, 0) == 0, "bip null");
    }

    /* Components + topo. */
    {
        int eu[] = {0, 2}, ev[] = {1, 3};
        ASSERT_TRUE(&ctx, ml_components(4, eu, ev, 2) == 2, "two components");
        ASSERT_TRUE(&ctx, ml_components(4, eu, ev, 0) == 4, "isolated");
        ASSERT_TRUE(&ctx, ml_components(4, NULL, NULL, 0) == -1, "comp null");
    }
    {
        int eu[] = {0, 0, 1}, ev[] = {1, 2, 2};
        int out[3] = {0, 0, 0};
        int c = ml_topo_sort(3, eu, ev, 3, out);
        ASSERT_TRUE(&ctx, c == 3, "topo count");
        ASSERT_TRUE(&ctx, out[0] == 0 && out[2] == 2, "topo order ends");
        {
            int cu[] = {0, 1}, cv[] = {1, 0};
            ASSERT_TRUE(&ctx, ml_topo_sort(2, cu, cv, 2, out) == -1, "topo cycle");
        }
    }

    /* Invalid + negative weights. */
    {
        int eu[] = {0}, ev[] = {1};
        double w[] = {1.0}, wn[] = {-1.0};
        double d[2], tot = 0;
        ASSERT_TRUE(&ctx, ml_dijkstra(2, eu, ev, wn, 1, 0, d) == ML_ERR_INVALID_ARG, "dij negative");
        ASSERT_TRUE(&ctx, ml_dijkstra(2, eu, ev, w, 1, 5, d) == ML_ERR_INVALID_ARG, "dij bad src");
        ASSERT_TRUE(&ctx, ml_kruskal_mst(2, eu, ev, wn, 1, &tot) == ML_ERR_INVALID_ARG, "kr negative");
    }

    return ml_test_summary(&ctx);
}
