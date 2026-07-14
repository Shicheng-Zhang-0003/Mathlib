#include <stdio.h>
#include "ml_core.h"
#include "ml_tensor.h"
#include "ml_linalg.h"
#include "ml_types.h"

int main() {
    printf("=== MathLib v11A1: Zero-Alloc Tensor Solver ===\n\n");

    // 1. Define a client-managed scratchpad (No malloc!)
    char scratchpad[4096];
    ml_workspace_t ws = { .storage = scratchpad, .size_bytes = sizeof(scratchpad), .used_bytes = 0 };
    ml_workspace_init(&ws);

    // 2. Define a simple 2x2 system: Ax = b
    // [ 2  0 ] [ x1 ]   [ 4 ]
    // [ 0  2 ] [ x2 ] = [ 6 ]
    double A_data[4] = {2.0, 0.0, 0.0, 2.0};
    double b_data[2] = {4.0, 6.0};
    double x_data[2] = {0.0, 0.0};

    ml_tensor_view_t A_view = ml_tensor_view(A_data, 2, 2);

    // 3. Solve using the zero-alloc engine
    ml_status_t status = ml_solve_v10(A_view, b_data, x_data, &ws);

    if (status == ML_SUCCESS) {
        printf("✅ Success! Solution:\n");
        printf("   x1 = %.2f\n", x_data[0]);
        printf("   x2 = %.2f\n", x_data[1]);
    } else {
        printf("❌ Solver failed with status: %d\n", status);
    }

    return 0;
}
