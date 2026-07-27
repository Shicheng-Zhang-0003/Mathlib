#include "test_harness.h"
#include "ml_tensor.h"
#include "ml_linalg.h"
#include "ml_types.h"

int main() {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Linear Algebra");
    printf("=== Linear Algebra Tests ===\n");

    char scratchpad[4096];
    ml_workspace_t ws;
    ml_workspace_init(&ws, scratchpad, sizeof(scratchpad));

    /* GLM Nitpick Fix: Use a non-trivial, non-symmetric matrix to actually
     * test Gaussian elimination and backward substitution, rather than a
     * trivial diagonal matrix that requires no pivoting. */
    double A[4] = {1, 2, 3, 4};
    double b[2] = {5, 11};
    double x[2] = {0};

    ml_tensor_view_t A_view = ml_tensor_view(A, 2, 2);
    int status = ml_solve(A_view, b, x, &ws);

    ASSERT_TRUE(&ctx, status == ML_SUCCESS, "Solver status");
    /* Solution to [1 2; 3 4] * [x1; x2] = [5; 11] is x1 = 1, x2 = 2 */
    ASSERT_NEAR(&ctx, x[0], 1.0, 1e-9, "x[0]");
    ASSERT_NEAR(&ctx, x[1], 2.0, 1e-9, "x[1]");

    return ml_test_summary(&ctx);
}
