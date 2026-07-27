/* v11S CLOSURE IP-20: edge statistics tests */
#include "test_harness.h"
#include "ml_statistics.h"

int main(void) {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Edge Statistics");

    double data[3] = {1.0, 2.0, 3.0};

    ASSERT_TRUE(&ctx, ml_isnan(ml_mean(NULL, 3)), "mean rejects NULL");
    ASSERT_TRUE(&ctx, ml_isnan(ml_mean(data, 0)), "mean rejects n <= 0");

    ASSERT_TRUE(&ctx, ml_isnan(ml_variance(NULL, 3)), "variance rejects NULL");
    ASSERT_TRUE(&ctx, ml_isnan(ml_variance(data, 0)), "variance rejects n <= 0");
    ASSERT_TRUE(&ctx, ml_isnan(ml_stddev(data, -1)), "stddev rejects n < 0");

    double m = 0.0;
    double b = 0.0;

    ml_linear_regression(NULL, data, 3, &m, &b);
    ASSERT_TRUE(&ctx, ml_isnan(m) && ml_isnan(b), "linear regression rejects NULL x");

    m = 0.0;
    b = 0.0;
    ml_linear_regression(data, data, 0, &m, &b);
    ASSERT_TRUE(&ctx, ml_isnan(m) && ml_isnan(b), "linear regression rejects n <= 0");

    ASSERT_TRUE(&ctx, ml_isnan(ml_binomial_pmf(-1, 0, 0.5)), "binomial rejects n < 0");
    ASSERT_TRUE(&ctx, ml_isnan(ml_binomial_pmf(5, 6, 0.5)), "binomial rejects k > n");
    ASSERT_TRUE(&ctx, ml_isnan(ml_binomial_pmf(5, 2, -0.1)), "binomial rejects p < 0");

    ASSERT_NEAR(&ctx, ml_binomial_pmf(5, 0, 0.0), 1.0, 1e-15, "binomial p=0,k=0");
    ASSERT_NEAR(&ctx, ml_binomial_pmf(5, 5, 1.0), 1.0, 1e-15, "binomial p=1,k=n");
    ASSERT_NEAR(&ctx, ml_binomial_pmf(5, 1, 0.0), 0.0, 1e-15, "binomial p=0,k>0");

    ASSERT_TRUE(&ctx, ml_isnan(ml_normal_pdf(0.0, 0.0, 0.0)), "normal_pdf rejects sigma <= 0");
    ASSERT_TRUE(&ctx, ml_normal_pdf(0.0, 0.0, ml_make_inf(0)) == 0.0, "normal_pdf infinite sigma -> 0");

    return ml_test_summary(&ctx);
}
