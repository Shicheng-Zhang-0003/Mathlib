/* v11S CLOSURE IP-20: edge linalg tests */
#include "test_harness.h"
#include "ml_tensor.h"
#include "ml_linalg.h"
#include "ml_types.h"

int main(void) {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Edge Linear Algebra");

    char scratch[4096];
    ml_workspace_t ws;
    ml_workspace_init(&ws, scratch, sizeof(scratch));

    ml_tensor_view_t Anull;
    Anull.data = NULL;
    Anull.rows = 2;
    Anull.cols = 2;

    double b[2] = {1.0, 1.0};
    double x[2] = {0.0, 0.0};

    ASSERT_TRUE(&ctx, ml_solve(Anull, b, x, &ws) == ML_ERR_INVALID_ARG, "solve rejects NULL A");

    double Adummy[1] = {0.0};
    ml_tensor_view_t Abad = ml_tensor_view(Adummy, 0, 0);
    ASSERT_TRUE(&ctx, ml_solve(Abad, b, x, &ws) == ML_ERR_INVALID_ARG, "solve rejects bad dimensions");

    char tiny_scratch[16];
    ml_workspace_t tiny_ws;
    ml_workspace_init(&tiny_ws, tiny_scratch, sizeof(tiny_scratch));

    double A2[4] = {2.0, 0.0, 0.0, 2.0};
    double b2[2] = {4.0, 6.0};
    double x2[2] = {0.0, 0.0};
    ml_tensor_view_t A2v = ml_tensor_view(A2, 2, 2);

    ASSERT_TRUE(&ctx, ml_solve(A2v, b2, x2, &tiny_ws) == ML_ERR_WORKSPACE, "solve detects workspace exhaustion");

    double Atiny[4] = {1e-300, 0.0, 0.0, 1e-300};
    double btiny[2] = {1e-300, 2e-300};
    double xtiny[2] = {0.0, 0.0};
    ml_tensor_view_t Atinyv = ml_tensor_view(Atiny, 2, 2);

    ml_workspace_reset(&ws);
    int status = ml_solve(Atinyv, btiny, xtiny, &ws);

    ASSERT_TRUE(&ctx, status == ML_SUCCESS, "tiny scaled matrix solve succeeds");
    ASSERT_NEAR(&ctx, xtiny[0], 1.0, 1e-12, "tiny scaled x[0]");
    ASSERT_NEAR(&ctx, xtiny[1], 2.0, 1e-12, "tiny scaled x[1]");

    double Az[4] = {0.0, 0.0, 0.0, 0.0};
    double bz[2] = {1.0, 1.0};
    double xz[2] = {0.0, 0.0};
    ml_tensor_view_t Azv = ml_tensor_view(Az, 2, 2);

    ml_workspace_reset(&ws);
    ASSERT_TRUE(&ctx, ml_solve(Azv, bz, xz, &ws) == ML_ERR_SINGULAR, "zero matrix is singular");

    ml_matvec(Anull, b, x);
    ASSERT_TRUE(&ctx, 1, "matvec NULL survived");

    return ml_test_summary(&ctx);
}
