#include "ml_compiler.h"
#include "ml_complex.h"
#include "internal/hypot.h"
#include "internal/pow_util.h"

/* v11S CLOSURE IP-4: overflow-safe complex geometry */

ML_API double ml_cplx_abs(cplx a) {
    return ml_hypot_internal(a.real, a.imag);
}

ML_API double ml_cplx_arg(cplx a) {
    /* MATHLIB_CLOSURE_P0_COMPLEX_ARG_ATAN2 */
    return ml_atan2(a.imag, a.real);
}

/* Smith's method for complex division */
ML_API cplx ml_cplx_div(cplx a, cplx b) {
    double ratio, denom;

    if (ml_fabs(b.real) >= ml_fabs(b.imag)) {
        if (b.real == 0.0) return (cplx){ml_make_nan(), ml_make_nan()};
        ratio = b.imag / b.real;
        denom = b.real + ratio * b.imag;
        if (denom == 0.0) return (cplx){ml_make_nan(), ml_make_nan()};
        return (cplx){
            (a.real + a.imag * ratio) / denom,
            (a.imag - a.real * ratio) / denom
        };
    } else {
        if (b.imag == 0.0) return (cplx){ml_make_nan(), ml_make_nan()};
        ratio = b.real / b.imag;
        denom = b.imag + ratio * b.real;
        if (denom == 0.0) return (cplx){ml_make_nan(), ml_make_nan()};
        return (cplx){
            (a.real * ratio + a.imag) / denom,
            (a.imag * ratio - a.real) / denom
        };
    }
}

ML_API cplx ml_cplx_exponential(cplx a) {
    if (ml_isnan(a.real) || ml_isnan(a.imag)) {
        return (cplx){ml_make_nan(), ml_make_nan()};
    }

    if (ml_isinf(a.imag)) {
        return (cplx){ml_make_nan(), ml_make_nan()};
    }

    if (ml_isinf(a.real)) {
        double c = ml_cos(a.imag);
        double s = ml_sin(a.imag);

        if (a.real < 0.0) {
            return (cplx){ml_copysign(0.0, c), ml_copysign(0.0, s)};
        }

        double re = (c == 0.0) ? ml_copysign(0.0, c)
                               : ml_copysign(ml_make_inf(0), c);
        double im = (s == 0.0) ? ml_copysign(0.0, s)
                               : ml_copysign(ml_make_inf(0), s);
        return (cplx){re, im};
    }

    double mag = ml_exp(a.real);

    if (ml_isinf(mag)) {
        double c = ml_cos(a.imag);
        double s = ml_sin(a.imag);

        double re = (c == 0.0) ? ml_copysign(0.0, c)
                               : ml_copysign(ml_make_inf(0), c);
        double im = (s == 0.0) ? ml_copysign(0.0, s)
                               : ml_copysign(ml_make_inf(0), s);
        return (cplx){re, im};
    }

    double c = ml_cos(a.imag);
    double s = ml_sin(a.imag);

    return (cplx){mag * c, mag * s};
}

ML_API cplx ml_cplx_logarithm(cplx a) {
    return (cplx){ml_log(ml_cplx_abs(a)), ml_cplx_arg(a)};
}

ML_API cplx ml_cplx_sqrt(cplx a) {
    /* Stable: r=|a|, re=sqrt((r+re)/2), im=sign*sqrt((r-re)/2). No overflow. */
    if (ml_isnan(a.real) || ml_isnan(a.imag)) {
        double n = ml_make_nan();
        return (cplx){n, n};
    }
    if (ml_isinf(a.real) || ml_isinf(a.imag)) {
        /* Infinite magnitude dominates: check imag-Inf first so
         * (-Inf+Infi) gives (Inf,Inf), not (0,Inf). */
        if (ml_isinf(a.imag)) {
            if (ml_isinf(a.real)) {
                if (a.real > 0.0) {
                    return (cplx){ml_make_inf(0), ml_copysign(ml_make_inf(0), a.imag)};
                }
                return (cplx){ml_make_inf(0), ml_copysign(ml_make_inf(0), a.imag)};
            }
            return (cplx){ml_make_inf(0), ml_copysign(ml_make_inf(0), a.real)};
        }
        if (ml_isinf(a.real)) {
            if (a.real > 0.0) {
                double m = ml_fabs(a.imag);
                (void)m;
                return (cplx){ml_make_inf(0), ml_copysign(ml_make_inf(0), a.imag)};
            }
            return (cplx){0.0, ml_copysign(ml_make_inf(0), a.imag)};
        }
        return (cplx){ml_make_inf(0), ml_copysign(ml_make_inf(0), a.real)};
    }
    {
        double r = ml_cplx_abs(a);
        if (r == 0.0) return (cplx){ml_copysign(0.0, a.real), ml_copysign(0.0, a.imag)};
        {
            double re = ml_sqrt((r + a.real) * 0.5);
            double im = ml_sqrt((r - a.real) * 0.5);
            if (ml_signbit(a.imag)) im = -im;
            return (cplx){re, im};
        }
    }
}

/* MATHLIB_CLOSURE_P0_COMPLEX_POWER_TREE */
static int ml_cplx_is_nan(cplx z) {
    return ml_isnan(z.real) || ml_isnan(z.imag);
}

static int ml_cplx_is_zero(cplx z) {
    return z.real == 0.0 && z.imag == 0.0;
}

static int ml_cplx_is_one(cplx z) {
    return z.real == 1.0 && z.imag == 0.0;
}

static int ml_cplx_is_pure_real_inf(cplx z) {
    return ml_isinf(z.real) && z.imag == 0.0;
}

static cplx ml_cplx_make_nan(void) {
    double n = ml_make_nan();
    return (cplx){n, n};
}

static cplx ml_cplx_make_zero(void) {
    return (cplx){0.0, 0.0};
}

static cplx ml_cplx_make_one(void) {
    return (cplx){1.0, 0.0};
}

static cplx ml_cplx_make_pos_inf(void) {
    return (cplx){ml_make_inf(0), 0.0};
}

static cplx ml_cplx_make_neg_inf(void) {
    return (cplx){-ml_make_inf(0), 0.0};
}

static cplx ml_cplx_make_neg_zero(void) {
    return (cplx){ml_copysign(0.0, -1.0), 0.0};
}

ML_API cplx ml_cplx_power(cplx a, cplx b) {
    int a_nan = ml_cplx_is_nan(a);
    int b_nan = ml_cplx_is_nan(b);

    int a_zero = ml_cplx_is_zero(a);
    int a_one = ml_cplx_is_one(a);
    int b_zero = ml_cplx_is_zero(b);

    int a_inf = ml_isinf(a.real) || ml_isinf(a.imag);
    int b_inf = ml_isinf(b.real) || ml_isinf(b.imag);

    if (b_zero) {
        return ml_cplx_make_one();
    }

    if (a_one) {
        return ml_cplx_make_one();
    }

    if (a_nan) {
        return ml_cplx_make_nan();
    }

    if (b_nan) {
        return ml_cplx_make_nan();
    }

    if (a_zero) {
        if (b_inf) {
            if (ml_isinf(b.imag)) {
                return ml_cplx_make_nan();
            }

            if (b.real > 0.0) {
                return ml_cplx_make_zero();
            }

            if (b.real < 0.0) {
                return ml_cplx_make_pos_inf();
            }

            return ml_cplx_make_nan();
        }

        if (b.real > 0.0) {
            return ml_cplx_make_zero();
        }

        if (b.real < 0.0) {
            return ml_cplx_make_pos_inf();
        }

        return ml_cplx_make_nan();
    }

    if (a_inf) {
        if (b_inf) {
            /* Pure-real infinite base and exponent: resolve by magnitude.
             * +Inf^+Inf=+Inf, +Inf^-Inf=+0; -Inf follows odd/even sign
             * when b is a finite odd integer, else magnitude rule. */
            if (ml_cplx_is_pure_real_inf(a) && b.imag == 0.0) {
                if (b.real > 0.0) {
                    if (a.real < 0.0 && ml_is_odd_integer_double(b.real)) {
                        return ml_cplx_make_neg_inf();
                    }
                    return ml_cplx_make_pos_inf();
                }
                if (b.real < 0.0) {
                    if (a.real < 0.0 && ml_is_odd_integer_double(b.real)) {
                        return ml_cplx_make_neg_zero();
                    }
                    return ml_cplx_make_zero();
                }
            }
            /* Non-pure-real infinities: fall through to exp(b*log a)
             * rather than unconditionally returning NaN. */
            if (!ml_cplx_is_pure_real_inf(a) || b.imag != 0.0) {
                cplx log_a = ml_cplx_logarithm(a);
                return ml_cplx_exponential((cplx){
                    b.real * log_a.real - b.imag * log_a.imag,
                    b.real * log_a.imag + b.imag * log_a.real
                });
            }
            return ml_cplx_make_nan();
        }

        if (!ml_cplx_is_pure_real_inf(a)) {
            return ml_cplx_make_nan();
        }

        if (b.imag != 0.0) {
            return ml_cplx_make_nan();
        }

        if (b.real > 0.0) {
            if (a.real < 0.0 && ml_is_odd_integer_double(b.real)) {
                return ml_cplx_make_neg_inf();
            }
            return ml_cplx_make_pos_inf();
        }

        if (b.real < 0.0) {
            if (a.real < 0.0 && ml_is_odd_integer_double(b.real)) {
                return ml_cplx_make_neg_zero();
            }
            return ml_cplx_make_zero();
        }

        return ml_cplx_make_one();
    }

    if (b_inf) {
        /* Finite nonzero a with infinite b: resolve via exp(b*log a)
         * instead of unconditional NaN. |a|>1 -> Inf, |a|<1 -> 0,
         * |a|==1 -> indeterminate (NaN). Zero-guard the complex product:
         * IEEE Inf*0=NaN would poison the real-real case
         * (e.g. cpow(2,Inf+0i) imag must be 0, not NaN). */
        cplx log_a = ml_cplx_logarithm(a);
        if (!ml_isfinite(log_a.real) || !ml_isfinite(log_a.imag)) {
            return ml_cplx_make_nan();
        }
        if (log_a.real == 0.0) {
            return ml_cplx_make_nan();
        }
        {
            double re, im;
            if (b.imag == 0.0 && log_a.imag == 0.0) {
                re = b.real * log_a.real;
                im = 0.0;
            } else if (b.imag == 0.0) {
                re = b.real * log_a.real;
                im = b.real * log_a.imag;
            } else if (log_a.imag == 0.0 && log_a.real != 0.0) {
                re = b.real * log_a.real;
                im = b.imag * log_a.real;
            } else {
                re = b.real * log_a.real - b.imag * log_a.imag;
                im = b.real * log_a.imag + b.imag * log_a.real;
            }
            return ml_cplx_exponential((cplx){re, im});
        }
    }

    cplx log_a = ml_cplx_logarithm(a);

    return ml_cplx_exponential((cplx){
        b.real * log_a.real - b.imag * log_a.imag,
        b.real * log_a.imag + b.imag * log_a.real
    });
}

ML_API cplx ml_cplx_nth_root(cplx a, int n, int k) {
    double nn = ml_make_nan();
    if (ML_UNLIKELY(n <= 0)) return (cplx){nn, nn};
    {
        int kk = ((k % n) + n) % n;
        double r = ml_cplx_abs(a);
        double th = ml_cplx_arg(a);
        if (!ml_isfinite(r) || ml_isnan(th)) {
            /* Infinite magnitude -> infinite root (angle moot).
             * Use signed Inf/zero to avoid Inf*0=NaN. */
            if (ml_isinf(r)) {
                double ph = 0.0;
                /* Preserve direction when defined. */
                if (ml_isfinite(th)) ph = (th + 2.0 * ML_PI * (double)kk) / (double)n;
                {
                    double c = ml_cos(ph), s = ml_sin(ph);
                    double re = (c == 0.0) ? ml_copysign(0.0, c) : ml_copysign(ml_make_inf(0), c);
                    double im = (s == 0.0) ? ml_copysign(0.0, s) : ml_copysign(ml_make_inf(0), s);
                    return (cplx){re, im};
                }
            }
            return (cplx){nn, nn};
        }
        if (r == 0.0) return (cplx){0.0, 0.0};
        {
            double m = ml_exp(ml_log(r) / (double)n);
            double ph = (th + 2.0 * ML_PI * (double)kk) / (double)n;
            return (cplx){m * ml_cos(ph), m * ml_sin(ph)};
        }
    }
}
