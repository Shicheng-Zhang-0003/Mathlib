#include <time.h>
#include "test_harness.h"

/* ============================================================================
 * v11S PUBLIC API SMOKE TEST
 *
 * PURGE NOTES:
 * - Eradicated 7 copy-pasted ml_exp/ml_log assertions.
 * - Removed dead code (test_ode_exp is now actually used).
 * - Removed embedded benchmark (had no assertions, inflated runtime).
 * - Removed internal/ header includes (tests public API routing instead).
 * - Transitioned to ml_test_ctx_t (no more global linker bombs).
 * - Updated workspace init to v11S ABI.
 * ========================================================================== */

#include "ml_combinatorics.h"
#include "ml_quadratics.h"
#include "ml_integral.h"
#include "ml_trig.h"
#include "ml_exp_log.h"
#include "ml_numerical.h"
#include "ml_polynomial.h"
#include "ml_complex.h"
#include "ml_statistics.h"
#include "ml_ode.h"
#include "ml_optimization.h"
#include "fft.h"
#include "bitwise_fp.h"
#include "fast_math.h"
#include "ml_tensor.h"
#include "ml_linalg.h"
#include "ml_fixed_point.h"
#include "profiles.h"
#include "ml_quaternion.h"
#include "simd.h"
#include "ieee754.h"

/* Test Helper Functions */
static double test_f_quad(double x) { return x * x - 4.0; }
static double test_df_quad(double x) { return 2.0 * x; }
static double test_ode_exp(double t, double y) { (void)t; return y; }
static double test_opt_parabola(double x) { return (x - 3.0) * (x - 3.0) + 1.0; }

int main() {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Monolithic Smoke Test");
    printf("=== MathLib v11S: Monolithic Smoke Test ===\n");

    /* 1. Combinatorics */
    printf("--- Combinatorics ---\n");
    ASSERT_TRUE(&ctx, ml_factorial(0) == 1, "0!");
    ASSERT_TRUE(&ctx, ml_factorial(5) == 120, "5!");
    ASSERT_TRUE(&ctx, ml_ncr(5, 2) == 10, "5 nCr 2");
    ASSERT_TRUE(&ctx, ml_npr(5, 2) == 20, "5 nPr 2");

    /* 2. Quadratics */
    printf("--- Quadratics ---\n");
    ASSERT_NEAR(&ctx, ml_formula_pos(1.0, 3.0, 2.0), -1.0, 1e-14, "formula_pos");
    ASSERT_NEAR(&ctx, ml_formula_neg(1.0, 3.0, 2.0), -2.0, 1e-14, "formula_neg");

    /* 3. Integral & Gamma */
    printf("--- Integral & Gamma ---\n");
    ASSERT_NEAR(&ctx, ml_integral_traditional(0, 1, 2.0, 0, 0.001), 1.0/3.0, 1e-3, "integral x^2");
    ASSERT_NEAR(&ctx, ml_gamma_new(1.0), 1.0, 1e-3, "gamma(1)");
    ASSERT_NEAR(&ctx, ml_gamma_new(4.0), 6.0, 1e-2, "gamma(4)");
    ASSERT_NEAR(&ctx, ml_gamma_new(5.0), 24.0, 5e-2, "gamma(5)");

    /* 4. Trigonometry */
    printf("--- Trigonometry ---\n");
    ASSERT_NEAR(&ctx, ml_sin(0.0), 0.0, 1e-14, "sin(0)");
    ASSERT_NEAR(&ctx, ml_cos(0.0), 1.0, 1e-14, "cos(0)");
    ASSERT_NEAR(&ctx, ml_tan(ML_PI / 4.0), 1.0, 1e-14, "tan(pi/4)");

    /* 5. Exponential & Log */
    printf("--- Exponential & Log ---\n");
    ASSERT_NEAR(&ctx, ml_exp(0.0), 1.0, 1e-15, "exp(0)");
    ASSERT_NEAR(&ctx, ml_log(1.0), 0.0, 1e-15, "log(1)");
    ASSERT_NEAR(&ctx, ml_pow(2.0, 3.0), 8.0, 1e-14, "pow(2,3)");

    /* 6. Numerical Methods */
    printf("--- Numerical Methods ---\n");
    ASSERT_NEAR(&ctx, ml_newton_raphson(test_f_quad, test_df_quad, 3.0, 1e-12, 100), 2.0, 1e-9, "newton_raphson");
    ASSERT_NEAR(&ctx, ml_bisection(test_f_quad, 1.0, 3.0, 1e-12, 100), 2.0, 1e-9, "bisection");

    /* 7. Polynomial */
    printf("--- Polynomial ---\n");
    double coeffs_quad[] = {-4.0, 0.0, 1.0};
    ASSERT_NEAR(&ctx, ml_polynomial_eval(coeffs_quad, 2, 2.0), 0.0, 1e-14, "poly_eval");

    /* 8. Complex */
    printf("--- Complex ---\n");
    cplx a = {1.0, 2.0};
    cplx b = {3.0, 4.0};
    cplx r = ml_cplx_add(a, b);
    ASSERT_NEAR(&ctx, r.real, 4.0, 1e-14, "cplx_add real");
    ASSERT_NEAR(&ctx, r.imag, 6.0, 1e-14, "cplx_add imag");

    /* 9. Statistics */
    printf("--- Statistics ---\n");
    double data[] = {2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0};
    ASSERT_NEAR(&ctx, ml_mean(data, 8), 5.0, 1e-14, "mean");
    ASSERT_NEAR(&ctx, ml_variance(data, 8), 4.0, 1e-14, "variance");

    /* 10. ODE (FIXED: Actually tests ODE now, instead of copy-pasted exp/log) */
    printf("--- ODE ---\n");
    ASSERT_NEAR(&ctx, ml_ode_euler(test_ode_exp, 0.0, 1.0, 0.001, 1000), 2.71692, 1e-2, "ode_euler");

    /* 11. Optimization */
    printf("--- Optimization ---\n");
    ASSERT_NEAR(&ctx, ml_optimize_golden(test_opt_parabola, -10.0, 10.0, 1e-5, 100), 3.0, 1e-4, "golden_section");

    /* 12. FFT */
    printf("--- FFT ---\n");
    cplx sig[4] = {{1.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};
    ml_fft_execute(sig, 4);
    ASSERT_NEAR(&ctx, sig[0].real, 1.0, 1e-9, "fft_dc");

    /* 13. Linear Algebra (FIXED: Uses v11S Workspace ABI) */
    printf("--- Linear Algebra ---\n");
    char scratchpad[4096];
    ml_workspace_t ws;
    ml_workspace_init(&ws, scratchpad, sizeof(scratchpad));

    double A[4] = {2.0, 0.0, 0.0, 2.0};
    double b_vec[2] = {4.0, 6.0};
    double x_vec[2] = {0.0, 0.0};
    ml_tensor_view_t A_view = ml_tensor_view(A, 2, 2);

    int status = ml_solve(A_view, b_vec, x_vec, &ws);
    ASSERT_TRUE(&ctx, status == ML_SUCCESS, "linalg_solve status");
    ASSERT_NEAR(&ctx, x_vec[0], 2.0, 1e-9, "linalg x[0]");

    /* 14. IEEE-754 Edge Cases */
    printf("--- IEEE-754 Edge Cases ---\n");
    ASSERT_TRUE(&ctx, ml_isnan(ml_asin(2.0)), "asin(2) is NaN");
    ASSERT_TRUE(&ctx, ml_isinf(ml_exp(1000.0)), "exp(1000) is Inf");

    return ml_test_summary(&ctx);
}
