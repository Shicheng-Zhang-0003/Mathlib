/* v11S AUDIT IP-3: complex special-case regression tests */
#include "test_harness.h"
#include "ml_complex.h"
#include "ml_core.h"

int main(void) {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "Edge Audit IP-3");

    cplx zero = {0.0, 0.0};
    cplx one = {1.0, 0.0};
    cplx mone = {-1.0, 0.0};
    cplx iunit = {0.0, 1.0};
    cplx niunit = {0.0, -1.0};
    cplx nanv = {ml_make_nan(), 0.0};
    cplx inf0 = {ml_make_inf(0), 0.0};
    cplx ninf0 = {-ml_make_inf(0), 0.0};
    cplx big = {1e308, 1e308};

    cplx p, e, l, d;

    /* abs / arg */
    ASSERT_TRUE(&ctx, ml_isfinite(ml_cplx_abs(big)), "abs overflow-safe");
    ASSERT_NEAR(&ctx, ml_cplx_abs(big), 1.4142135623730951e308, 1e296, "abs value");
    ASSERT_NEAR(&ctx, ml_cplx_abs(zero), 0.0, 1e-15, "abs zero");
    ASSERT_TRUE(&ctx, ml_isnan(ml_cplx_abs(nanv)), "abs NaN");
    ASSERT_TRUE(&ctx, ml_isinf(ml_cplx_abs(inf0)), "abs inf");

    ASSERT_NEAR(&ctx, ml_cplx_arg(one), 0.0, 1e-15, "arg(1)");
    ASSERT_NEAR(&ctx, ml_cplx_arg(iunit), ML_PI / 2.0, 1e-15, "arg(i)");
    ASSERT_NEAR(&ctx, ml_cplx_arg(mone), ML_PI, 1e-15, "arg(-1)");
    ASSERT_NEAR(&ctx, ml_cplx_arg(niunit), -ML_PI / 2.0, 1e-15, "arg(-i)");
    ASSERT_NEAR(&ctx, ml_cplx_arg(zero), 0.0, 1e-15, "arg(0)");
    ASSERT_TRUE(&ctx, ml_isnan(ml_cplx_arg(nanv)), "arg(NaN)");

    /* div */
    d = ml_cplx_div((cplx){1.0, 2.0}, (cplx){3.0, 4.0});
    ASSERT_NEAR(&ctx, d.real, 0.44, 1e-12, "div real");
    ASSERT_NEAR(&ctx, d.imag, 0.08, 1e-12, "div imag");

    d = ml_cplx_div(one, iunit);
    ASSERT_NEAR(&ctx, d.real, 0.0, 1e-12, "1/i real");
    ASSERT_NEAR(&ctx, d.imag, -1.0, 1e-12, "1/i imag");

    d = ml_cplx_div(iunit, one);
    ASSERT_NEAR(&ctx, d.real, 0.0, 1e-12, "i/1 real");
    ASSERT_NEAR(&ctx, d.imag, 1.0, 1e-12, "i/1 imag");

    d = ml_cplx_div(mone, one);
    ASSERT_NEAR(&ctx, d.real, -1.0, 1e-12, "-1/1 real");
    ASSERT_NEAR(&ctx, d.imag, 0.0, 1e-12, "-1/1 imag");

    d = ml_cplx_div(one, mone);
    ASSERT_NEAR(&ctx, d.real, -1.0, 1e-12, "1/-1 real");
    ASSERT_NEAR(&ctx, d.imag, 0.0, 1e-12, "1/-1 imag");

    d = ml_cplx_div(zero, one);
    ASSERT_NEAR(&ctx, d.real, 0.0, 1e-12, "0/1 real");
    ASSERT_NEAR(&ctx, d.imag, 0.0, 1e-12, "0/1 imag");

    ASSERT_TRUE(&ctx, ml_isnan(ml_cplx_div(one, zero).real) &&
                      ml_isnan(ml_cplx_div(one, zero).imag), "div by zero NaN");
    ASSERT_TRUE(&ctx, ml_isnan(ml_cplx_div(zero, zero).real) &&
                      ml_isnan(ml_cplx_div(zero, zero).imag), "0/0 NaN");
    ASSERT_TRUE(&ctx, ml_isnan(ml_cplx_div(nanv, one).real) &&
                      ml_isnan(ml_cplx_div(nanv, one).imag), "NaN/1 NaN");
    ASSERT_TRUE(&ctx, ml_isnan(ml_cplx_div(one, nanv).real) &&
                      ml_isnan(ml_cplx_div(one, nanv).imag), "1/NaN NaN");

    /* exp */
    e = ml_cplx_exponential(zero);
    ASSERT_NEAR(&ctx, e.real, 1.0, 1e-12, "exp(0) real");
    ASSERT_NEAR(&ctx, e.imag, 0.0, 1e-12, "exp(0) imag");

    e = ml_cplx_exponential(one);
    ASSERT_NEAR(&ctx, e.real, ML_E, 1e-12, "exp(1) real");
    ASSERT_NEAR(&ctx, e.imag, 0.0, 1e-12, "exp(1) imag");

    e = ml_cplx_exponential((cplx){0.0, ML_PI});
    ASSERT_NEAR(&ctx, e.real, -1.0, 1e-12, "exp(i*pi) real");
    ASSERT_NEAR(&ctx, e.imag, 0.0, 1e-12, "exp(i*pi) imag");

    e = ml_cplx_exponential((cplx){1000.0, 0.0});
    ASSERT_TRUE(&ctx, ml_isinf(e.real) && e.real > 0.0, "exp(1000) inf");

    e = ml_cplx_exponential((cplx){-1000.0, 0.0});
    ASSERT_TRUE(&ctx, e.real == 0.0 && e.imag == 0.0, "exp(-1000) zero");

    e = ml_cplx_exponential(inf0);
    ASSERT_TRUE(&ctx, ml_isinf(e.real) && e.real > 0.0 && e.imag == 0.0, "exp(inf+0i)");

    e = ml_cplx_exponential(ninf0);
    ASSERT_TRUE(&ctx, e.real == 0.0 && e.imag == 0.0, "exp(-inf+0i)");

    e = ml_cplx_exponential(nanv);
    ASSERT_TRUE(&ctx, ml_isnan(e.real) && ml_isnan(e.imag), "exp(NaN)");

    e = ml_cplx_exponential((cplx){0.0, ml_make_inf(0)});
    ASSERT_TRUE(&ctx, ml_isnan(e.real), "exp(0+i*inf) NaN");

    /* log */
    l = ml_cplx_logarithm(one);
    ASSERT_NEAR(&ctx, l.real, 0.0, 1e-12, "log(1) real");
    ASSERT_NEAR(&ctx, l.imag, 0.0, 1e-12, "log(1) imag");

    l = ml_cplx_logarithm((cplx){ML_E, 0.0});
    ASSERT_NEAR(&ctx, l.real, 1.0, 1e-12, "log(e) real");
    ASSERT_NEAR(&ctx, l.imag, 0.0, 1e-12, "log(e) imag");

    l = ml_cplx_logarithm(zero);
    ASSERT_TRUE(&ctx, ml_isinf(l.real) && l.real < 0.0, "log(0) -inf");

    l = ml_cplx_logarithm(nanv);
    ASSERT_TRUE(&ctx, ml_isnan(l.real) && ml_isnan(l.imag), "log(NaN)");

    /* power: zero / one / NaN special cases */
    p = ml_cplx_power(zero, zero);
    ASSERT_TRUE(&ctx, p.real == 1.0 && p.imag == 0.0, "0^0 = 1");

    p = ml_cplx_power(one, nanv);
    ASSERT_TRUE(&ctx, p.real == 1.0 && p.imag == 0.0, "1^NaN = 1");

    p = ml_cplx_power(nanv, zero);
    ASSERT_TRUE(&ctx, p.real == 1.0 && p.imag == 0.0, "NaN^0 = 1");

    p = ml_cplx_power(nanv, one);
    ASSERT_TRUE(&ctx, ml_isnan(p.real) && ml_isnan(p.imag), "NaN^1 NaN");

    p = ml_cplx_power(zero, nanv);
    ASSERT_TRUE(&ctx, ml_isnan(p.real) && ml_isnan(p.imag), "0^NaN NaN");

    p = ml_cplx_power(one, one);
    ASSERT_TRUE(&ctx, p.real == 1.0 && p.imag == 0.0, "1^1 = 1");

    p = ml_cplx_power(one, zero);
    ASSERT_TRUE(&ctx, p.real == 1.0 && p.imag == 0.0, "1^0 = 1");

    p = ml_cplx_power(zero, one);
    ASSERT_TRUE(&ctx, p.real == 0.0 && p.imag == 0.0, "0^1 = 0");

    p = ml_cplx_power(zero, (cplx){-1.0, 0.0});
    ASSERT_TRUE(&ctx, ml_isinf(p.real) && p.real > 0.0, "0^-1 = inf");

    p = ml_cplx_power(zero, (cplx){1.0, 1.0});
    ASSERT_TRUE(&ctx, p.real == 0.0 && p.imag == 0.0, "0^(1+i) = 0");

    p = ml_cplx_power(zero, (cplx){-1.0, 1.0});
    ASSERT_TRUE(&ctx, ml_isinf(p.real) && p.real > 0.0, "0^(-1+i) = inf");

    p = ml_cplx_power(zero, (cplx){0.0, 1.0});
    ASSERT_TRUE(&ctx, ml_isnan(p.real) && ml_isnan(p.imag), "0^i NaN");

    p = ml_cplx_power(zero, (cplx){0.0, -1.0});
    ASSERT_TRUE(&ctx, ml_isnan(p.real) && ml_isnan(p.imag), "0^-i NaN");

    p = ml_cplx_power(zero, inf0);
    ASSERT_TRUE(&ctx, p.real == 0.0 && p.imag == 0.0, "0^inf = 0");

    p = ml_cplx_power(zero, ninf0);
    ASSERT_TRUE(&ctx, ml_isinf(p.real) && p.real > 0.0, "0^-inf = inf");

    p = ml_cplx_power(zero, (cplx){ml_make_inf(0), 1.0});
    ASSERT_TRUE(&ctx, p.real == 0.0 && p.imag == 0.0, "0^(inf+i) = 0");

    p = ml_cplx_power(zero, (cplx){-ml_make_inf(0), 1.0});
    ASSERT_TRUE(&ctx, ml_isinf(p.real) && p.real > 0.0, "0^(-inf+i) = inf");

    p = ml_cplx_power(zero, (cplx){0.0, ml_make_inf(0)});
    ASSERT_TRUE(&ctx, ml_isnan(p.real) && ml_isnan(p.imag), "0^(i*inf) NaN");

    p = ml_cplx_power(one, inf0);
    ASSERT_TRUE(&ctx, p.real == 1.0 && p.imag == 0.0, "1^inf = 1");

    p = ml_cplx_power(one, ninf0);
    ASSERT_TRUE(&ctx, p.real == 1.0 && p.imag == 0.0, "1^-inf = 1");

    p = ml_cplx_power(one, (cplx){0.0, ml_make_inf(0)});
    ASSERT_TRUE(&ctx, p.real == 1.0 && p.imag == 0.0, "1^(i*inf) = 1");

    p = ml_cplx_power(one, (cplx){ml_make_nan(), 0.0});
    ASSERT_TRUE(&ctx, p.real == 1.0 && p.imag == 0.0, "1^(NaN+0i) = 1");

    p = ml_cplx_power(one, (cplx){0.0, ml_make_nan()});
    ASSERT_TRUE(&ctx, p.real == 1.0 && p.imag == 0.0, "1^(0+iNaN) = 1");

    p = ml_cplx_power((cplx){1.0, ml_make_nan()}, zero);
    ASSERT_TRUE(&ctx, p.real == 1.0 && p.imag == 0.0, "(1+iNaN)^0 = 1");

    p = ml_cplx_power((cplx){1.0, ml_make_nan()}, one);
    ASSERT_TRUE(&ctx, ml_isnan(p.real) && ml_isnan(p.imag), "(1+iNaN)^1 NaN");

    p = ml_cplx_power((cplx){ml_make_nan(), 0.0}, zero);
    ASSERT_TRUE(&ctx, p.real == 1.0 && p.imag == 0.0, "(NaN+0i)^0 = 1");

    p = ml_cplx_power((cplx){ml_make_nan(), 0.0}, one);
    ASSERT_TRUE(&ctx, ml_isnan(p.real) && ml_isnan(p.imag), "(NaN+0i)^1 NaN");

    p = ml_cplx_power((cplx){0.0, ml_make_nan()}, zero);
    ASSERT_TRUE(&ctx, p.real == 1.0 && p.imag == 0.0, "(0+iNaN)^0 = 1");

    p = ml_cplx_power((cplx){0.0, ml_make_nan()}, one);
    ASSERT_TRUE(&ctx, ml_isnan(p.real) && ml_isnan(p.imag), "(0+iNaN)^1 NaN");

    p = ml_cplx_power(inf0, zero);
    ASSERT_TRUE(&ctx, p.real == 1.0 && p.imag == 0.0, "inf^0 = 1");

    p = ml_cplx_power(inf0, one);
    ASSERT_TRUE(&ctx, ml_isinf(p.real) && p.real > 0.0, "inf^1 = inf");

    p = ml_cplx_power(inf0, (cplx){-1.0, 0.0});
    ASSERT_TRUE(&ctx, p.real == 0.0 && p.imag == 0.0, "inf^-1 = 0");

    p = ml_cplx_power(inf0, iunit);
    ASSERT_TRUE(&ctx, ml_isnan(p.real) && ml_isnan(p.imag), "inf^i NaN");

    p = ml_cplx_power(inf0, (cplx){1.0, 1.0});
    ASSERT_TRUE(&ctx, ml_isnan(p.real) && ml_isnan(p.imag), "inf^(1+i) NaN");

    p = ml_cplx_power(inf0, (cplx){-1.0, 1.0});
    ASSERT_TRUE(&ctx, ml_isnan(p.real) && ml_isnan(p.imag), "inf^(-1+i) NaN");

    /* negative real powers through principal branch */
    p = ml_cplx_power(mone, one);
    ASSERT_NEAR(&ctx, p.real, -1.0, 1e-12, "(-1)^1 real");
    ASSERT_NEAR(&ctx, p.imag, 0.0, 1e-12, "(-1)^1 imag");

    p = ml_cplx_power(mone, (cplx){2.0, 0.0});
    ASSERT_NEAR(&ctx, p.real, 1.0, 1e-12, "(-1)^2 real");
    ASSERT_NEAR(&ctx, p.imag, 0.0, 1e-12, "(-1)^2 imag");

    p = ml_cplx_power(mone, (cplx){0.5, 0.0});
    ASSERT_NEAR(&ctx, p.real, 0.0, 1e-12, "sqrt(-1) real");
    ASSERT_NEAR(&ctx, p.imag, 1.0, 1e-12, "sqrt(-1) imag");

    p = ml_cplx_power(mone, zero);
    ASSERT_TRUE(&ctx, p.real == 1.0 && p.imag == 0.0, "(-1)^0 = 1");

    p = ml_cplx_power(mone, nanv);
    ASSERT_TRUE(&ctx, ml_isnan(p.real) && ml_isnan(p.imag), "(-1)^NaN NaN");

    return ml_test_summary(&ctx);
}
