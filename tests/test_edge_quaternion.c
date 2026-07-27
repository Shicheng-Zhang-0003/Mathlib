/* v11S CLOSURE IP-20: edge quaternion tests */
#include "test_harness.h"
#include "ml_quaternion.h"

static double qnorm2(ml_quat q) {
    return q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z;
}

int main(void) {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Edge Quaternion");

    ml_quat zero = {0.0, 0.0, 0.0, 0.0};
    ml_quat nz = ml_quat_normalize(zero);
    ASSERT_TRUE(&ctx, ml_isnan(nz.w), "normalize zero quaternion gives NaN");

    ml_quat unit = {1.0, 0.0, 0.0, 0.0};
    ml_quat nunit = ml_quat_normalize(unit);
    ASSERT_NEAR(&ctx, nunit.w, 1.0, 1e-15, "normalize unit quaternion");

    ml_quat a = {1.0, 2.0, 3.0, 4.0};
    ml_quat b = {5.0, 6.0, 7.0, 8.0};
    ml_quat ab = ml_quat_mul(a, b);

    double na = qnorm2(a);
    double nb = qnorm2(b);
    double nab = qnorm2(ab);

    ASSERT_NEAR(&ctx, nab, na * nb, ml_fabs(na * nb) * 1e-12 + 1e-9, "norm(q1*q2) == norm(q1)*norm(q2)");

    ml_quat s0 = ml_quat_slerp(unit, (ml_quat){0.0, 1.0, 0.0, 0.0}, 0.0);
    ASSERT_NEAR(&ctx, s0.w, 1.0, 1e-12, "slerp t=0");

    ml_quat s1 = ml_quat_slerp(unit, (ml_quat){0.0, 1.0, 0.0, 0.0}, 1.0);
    ASSERT_NEAR(&ctx, s1.x, 1.0, 1e-12, "slerp t=1");

    return ml_test_summary(&ctx);
}
