#include "test_harness.h"
#include "tensor.h"
#include "linalg_v10.h"

int main() {
    printf("=== Linear Algebra Tests ===\n");
    char scratchpad[4096];
    ml_workspace_t ws = { scratchpad, sizeof(scratchpad), 0 };
    double A[4] = {2, 0, 0, 2};
    double b[2] = {4, 6};
    double x[2] = {0};
    ml_tensor_view_t A_view = ml_tensor_view(A, 2, 2);
    int status = ml_solve_v10(A_view, b, x, &ws);
    ASSERT_TRUE(status == 0, "Solver status");
    ASSERT_NEAR(x[0], 2.0, 1e-9, "x[0]");
    ASSERT_NEAR(x[1], 3.0, 1e-9, "x[1]");
    return test_harness_summary("Linalg");
}