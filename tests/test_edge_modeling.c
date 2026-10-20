/* Modeling pins: Markov, SIR, Leslie, M/M/1, conjugate Bayes, Lyapunov.
 * Ground truth: closed forms (4/7, e^-1, R0 formulas) + mpmath spot. */
#include "test_harness.h"
#include "ml_modeling.h"
#include "ml_core.h"
int main (void) {
    ml_test_ctx_t ctx;
    ml_test_init (&ctx, "Edge Modeling");

    /* Two-state chain [[.7,.3],[.4,.6]]: stationary (4/7, 3/7). */
    {
        double P [] = {0.7, 0.3, 0.4, 0.6};
        double pi [2] = {0, 0}, o [2] = {0, 0};
        ASSERT_TRUE (&ctx, ml_markov_is_stochastic (P, 2, 1e-12) == 1, "stochastic");
        ASSERT_TRUE (&ctx, ml_markov_is_stochastic ((double []) {0.5, 0.6, 0.5, 0.4}, 2, 1e-12) == 0, "not stochastic");
        ASSERT_TRUE (&ctx, ml_markov_step (P, (double []) {1, 0}, o, 2) == ML_SUCCESS, "step st");
        ASSERT_NEAR (&ctx, o [0], 0.7, 1e-12, "step o0");
        ASSERT_NEAR (&ctx, o [1], 0.3, 1e-12, "step o1");
        ASSERT_TRUE (&ctx, ml_markov_stationary (P, pi, 2, 10000, 1e-12) == ML_SUCCESS, "stat st");
        ASSERT_NEAR (&ctx, pi [0], 0.5714285714285714, 1e-9, "stat 4/7");
        ASSERT_NEAR (&ctx, pi [1], 0.4285714285714286, 1e-9, "stat 3/7");
    }

    /* SIR. */
    ASSERT_NEAR (&ctx, ml_sir_r0 (0.3, 0.1), 3.0, 1e-12, "R0");
    ASSERT_TRUE (&ctx, ml_isnan (ml_sir_r0 (0.3, 0.0)), "R0 zero gamma NaN");
    ASSERT_NEAR (&ctx, ml_sir_final_size (3.0), 0.9404797907, 1e-7, "Z R0=3");
    ASSERT_NEAR (&ctx, ml_sir_final_size (1.5), 0.5828116439, 1e-7, "Z R0=1.5");
    ASSERT_NEAR (&ctx, ml_sir_final_size (1.0), 0.0, 0.0, "Z R0=1");
    ASSERT_NEAR (&ctx, ml_sir_peak (3.0), 0.3004625704, 1e-7, "peak R0=3");
    ASSERT_NEAR (&ctx, ml_sir_peak (0.5), 0.0, 0.0, "peak subcritical");

    /* Leslie 2-class [[0.5,2],[0.4,0]]: lambda = (0.5+sqrt(3.45))/2. */
    {
        double L [] = {0.5, 2.0, 0.4, 0.0};
        double pop [2] = {0, 0}, lam = 0;
        ASSERT_TRUE (&ctx, ml_leslie_step (L, (double []) {10, 5}, pop, 2) == ML_SUCCESS, "leslie st");
        ASSERT_NEAR (&ctx, pop [0], 15.0, 1e-12, "leslie newborns");
        ASSERT_NEAR (&ctx, pop [1], 4.0, 1e-12, "leslie survivors");
        ASSERT_TRUE (&ctx, ml_leslie_growth (L, 2, &lam) == ML_SUCCESS, "growth st");
        ASSERT_NEAR (&ctx, lam, 1.1787090317, 1e-6, "leslie lambda");
    }
    /* Imprimitive [[0,2.5],[0.4,0]] oscillates: honest CONVERGENCE. */
    {
        double L [] = {0.0, 2.5, 0.4, 0.0};
        double lam = 0;
        ASSERT_TRUE (&ctx, ml_leslie_growth (L, 2, &lam) == ML_ERR_CONVERGENCE, "imprimitive converges-loud");
    }

    /* M/M/1: lambda=2, mu=3 -> rho=2/3, L=2, W=1, P(wait>1)=e^-1. */
    ASSERT_NEAR (&ctx, ml_mm1_rho (2, 3), 0.6666666666666666, 1e-12, "rho");
    ASSERT_NEAR (&ctx, ml_mm1_L (2, 3), 2.0, 1e-12, "L");
    ASSERT_NEAR (&ctx, ml_mm1_W (2, 3), 1.0, 1e-12, "W");
    ASSERT_NEAR (&ctx, ml_mm1_wait (2, 3, 1.0), 0.36787944117144233, 1e-12, "wait tail");
    ASSERT_TRUE (&ctx, ml_isnan (ml_mm1_L (3, 3)), "unstable NaN");
    ASSERT_TRUE (&ctx, ml_isnan (ml_mm1_W (-1, 3)), "negative NaN");

    /* Conjugate Bayes. */
    {
        double a1 = 0, b1 = 0, lo = 0, hi = 0;
        ASSERT_TRUE (&ctx, ml_beta_binomial_post (2, 2, 7, 10, &a1, &b1) == ML_SUCCESS, "bb st");
        ASSERT_NEAR (&ctx, a1, 9.0, 0.0, "bb a1");
        ASSERT_NEAR (&ctx, b1, 5.0, 0.0, "bb b1");
        ASSERT_NEAR (&ctx, ml_beta_post_mean (a1, b1), 0.6428571428571429, 1e-12, "post mean");
        ASSERT_TRUE (&ctx, ml_beta_credible (a1, b1, 0.05, &lo, &hi) == ML_SUCCESS, "cred st");
        ASSERT_TRUE (&ctx, lo > 0.35 && lo < 0.42 && hi > 0.83 && hi < 0.89, "cred interval sane");
        ASSERT_TRUE (&ctx, ml_beta_binomial_post (2, 2, 11, 10, &a1, &b1) == ML_ERR_INVALID_ARG, "succ>trials");
        ASSERT_TRUE (&ctx, ml_gamma_poisson_post (2, 1, 5, 2.0, &a1, &b1) == ML_SUCCESS, "gp st");
        ASSERT_NEAR (&ctx, a1, 7.0, 0.0, "gp a1");
        ASSERT_NEAR (&ctx, b1, 3.0, 0.0, "gp b1");
    }

    /* Chaos: r=2.5 -> -ln2 (fixed point), r=4 -> ln2 (fully chaotic). */
    ASSERT_NEAR (&ctx, ml_logistic_lyapunov (2.5, 0.3, 1000, 5000), -0.6931471805599453, 1e-3, "lyap stable");
    ASSERT_NEAR (&ctx, ml_logistic_lyapunov (4.0, 0.3, 1000, 5000), 0.6931471805599453, 1e-3, "lyap chaotic");
    ASSERT_TRUE (&ctx, ml_isnan (ml_logistic_lyapunov (5.0, 0.3, 100, 100)), "lyap r>4 NaN");

    return ml_test_summary (&ctx);
}
