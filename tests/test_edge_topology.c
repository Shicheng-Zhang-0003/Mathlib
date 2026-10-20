/* Doctorate topology pins: Euler, Betti F2, chain check, Bareiss. */
#include "test_harness.h"
#include "ml_topology.h"
#include "ml_core.h"
int main (void) {
    ml_test_ctx_t ctx;
    ml_test_init (&ctx, "Edge Topology");

    ASSERT_TRUE (&ctx, ml_euler_char ((int []) {4, 6, 4}, 2) == 2, "euler tet surface");
    ASSERT_TRUE (&ctx, ml_euler_char ((int []) {3, 3}, 1) == 0, "euler triangle loop");
    ASSERT_TRUE (&ctx, ml_euler_char (NULL, 2) == INT64_MIN, "euler null");

    /* Filled triangle: (1,0,0). */
    {
        int v [] = {0, 1, 2, 0, 1, 0, 2, 1, 2, 0, 1, 2};
        int off [] = {0, 1, 2, 3, 5, 7, 9, 12};
        int dim [] = {0, 0, 0, 1, 1, 1, 2};
        int b [3] = {0, 0, 0};
        ml_status_t st = ml_simplicial_betti_f2 (3, v, off, dim, 7, b, 2);
        ASSERT_TRUE (&ctx, st == ML_SUCCESS, "filled status");
        ASSERT_TRUE (&ctx, b [0] == 1 && b [1] == 0 && b [2] == 0, "filled betti");
    }
    /* Hollow triangle: (1,1,0). */
    {
        int v [] = {0, 1, 2, 0, 1, 0, 2, 1, 2};
        int off [] = {0, 1, 2, 3, 5, 7, 9};
        int dim [] = {0, 0, 0, 1, 1, 1};
        int b [3] = {0, 0, 0};
        ml_status_t st = ml_simplicial_betti_f2 (3, v, off, dim, 6, b, 2);
        ASSERT_TRUE (&ctx, st == ML_SUCCESS, "hollow status");
        ASSERT_TRUE (&ctx, b [0] == 1 && b [1] == 1 && b [2] == 0, "hollow betti");
    }
    /* Tetrahedron surface: (1,0,1). */
    {
        int v [] = {0, 1, 2, 3, 0, 1, 0, 2, 0, 3, 1, 2, 1, 3, 2, 3,
                    0, 1, 2, 0, 1, 3, 0, 2, 3, 1, 2, 3};
        int off [] = {0, 1, 2, 3, 4, 6, 8, 10, 12, 14, 16, 19, 22, 25, 28};
        int dim [] = {0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2};
        int b [4] = {0, 0, 0, 0};
        ml_status_t st = ml_simplicial_betti_f2 (4, v, off, dim, 14, b, 3);
        ASSERT_TRUE (&ctx, st == ML_SUCCESS, "tet status");
        ASSERT_TRUE (&ctx, b [0] == 1 && b [1] == 0 && b [2] == 1, "tet betti");
    }
    /* Non-closed: edge [0,1] without vertex [1]. */
    {
        int v [] = {0, 0, 1};
        int off [] = {0, 1, 3};
        int dim [] = {0, 1};
        int b [2] = {0, 0};
        ASSERT_TRUE (&ctx, ml_simplicial_betti_f2 (2, v, off, dim, 2, b, 1) == ML_ERR_INVALID_ARG, "non-closed");
    }

    /* Chain 0 -> C2 --[1 1]--> C1 --[1 1]^T--> C0 -> 0 is a complex
     * ([1 1] dot [1 1]^T = 0 mod 2) and acyclic: betti (0,0,0). */
    {
        int mats [] = {1, 1, 1, 1};
        int rows [] = {1, 2};
        int cols [] = {2, 1};
        int offm [] = {0, 2};
        int b [3] = {0, 0, 0};
        ASSERT_TRUE (&ctx, ml_chain_f2_is_complex (mats, rows, cols, offm, 2) == 1, "chain complex");
        ASSERT_TRUE (&ctx, ml_chain_betti_f2 (mats, rows, cols, offm, 2, b) == ML_SUCCESS, "chain betti st");
        ASSERT_TRUE (&ctx, b [0] == 0 && b [1] == 0 && b [2] == 0, "chain betti acyclic");
        {
            int bad [] = {1, 0, 1, 0};
            ASSERT_TRUE (&ctx, ml_chain_f2_is_complex (bad, rows, cols, offm, 2) == 0, "not a complex");
        }
    }
    /* Single fold map [1 1]: H = (0, 1). */
    {
        int mats [] = {1, 1};
        int rows [] = {1};
        int cols [] = {2};
        int offm [] = {0};
        int b [2] = {0, 0};
        ASSERT_TRUE (&ctx, ml_chain_betti_f2 (mats, rows, cols, offm, 1, b) == ML_SUCCESS, "fold betti st");
        ASSERT_TRUE (&ctx, b [0] == 0 && b [1] == 1, "fold betti");
    }

    /* Bareiss exactness. */
    {
        long long d = 0;
        ASSERT_TRUE (&ctx, ml_bareiss_det ((long long []) {1, 2, 3, 4}, 2, &d) == ML_SUCCESS && d == -2, "det -2");
        ASSERT_TRUE (&ctx, ml_bareiss_det ((long long []) {2, -1, 0, -1, 2, -1, 0, -1, 2}, 3, &d) == ML_SUCCESS && d == 4, "det 4");
        ASSERT_TRUE (&ctx, ml_bareiss_rank ((long long []) {1, 2, 2, 4}, 2, 2) == 1, "rank 1");
        ASSERT_TRUE (&ctx, ml_bareiss_rank ((long long []) {1, 0, 0, 1}, 2, 2) == 2, "rank 2");
        ASSERT_TRUE (&ctx, ml_bareiss_rank (NULL, 2, 2) == -1, "rank null");
    }

    return ml_test_summary (&ctx);
}
