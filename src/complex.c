#include "ml_compiler.h"
#include "ml_complex.h"
#include "internal/hypot.h"
#include "internal/pow_util.h"
/* v11S CLOSURE IP-4: overflow-safe complex geometry */

ML_API double ml_cplx_abs (cplx a) {
    return ml_hypot_internal (a.real, a.imag);
}
ML_API double ml_cplx_arg (cplx a) {
    /* MATHLIB_CLOSURE_P0_COMPLEX_ARG_ATAN2 */
    return ml_atan2 (a.imag, a.real);
}
/* Smith's method for complex division */
ML_API cplx ml_cplx_div (cplx a, cplx b) {
    double ratio, denom;

    if (ml_fabs (b.real) >= ml_fabs (b.imag)) {
        if (b.real == 0.0) return (cplx) {ml_make_nan (), ml_make_nan ()};
        ratio = b.imag / b.real;
        denom = b.real + ratio * b.imag;
        if (denom == 0.0) return (cplx) {ml_make_nan (), ml_make_nan ()};
        return (cplx) {
            (a.real + a.imag * ratio) / denom,
            (a.imag - a.real * ratio) / denom};
    } else {
        if (b.imag == 0.0) return (cplx) {ml_make_nan (), ml_make_nan ()};
        ratio = b.real / b.imag;
        denom = b.imag + ratio * b.real;
        if (denom == 0.0) return (cplx) {ml_make_nan (), ml_make_nan ()};
        return (cplx) {
            (a.real * ratio + a.imag) / denom,
            (a.imag * ratio - a.real) / denom};
    }
}
ML_API cplx ml_cplx_exponential (cplx a) {
    if (ml_isnan (a.real) || ml_isnan (a.imag)) {
        return (cplx) {ml_make_nan (), ml_make_nan ()};
    }

    if (ml_isinf (a.imag)) {
        return (cplx) {ml_make_nan (), ml_make_nan ()};
    }

    if (ml_isinf (a.real)) {
        double c = ml_cos (a.imag);
        double s = ml_sin (a.imag);

        if (a.real < 0.0) {
            return (cplx) {ml_copysign (0.0, c), ml_copysign (0.0, s)};
        }

        double re = (c == 0.0) ? ml_copysign (0.0, c)
                               : ml_copysign (ml_make_inf (0), c);
        double im = (s == 0.0) ? ml_copysign (0.0, s)
                               : ml_copysign (ml_make_inf (0), s);
        return (cplx) {re, im};
    }

    double mag = ml_exp (a.real);

    if (ml_isinf (mag)) {
        double c = ml_cos (a.imag);
        double s = ml_sin (a.imag);

        double re = (c == 0.0) ? ml_copysign (0.0, c)
                               : ml_copysign (ml_make_inf (0), c);
        double im = (s == 0.0) ? ml_copysign (0.0, s)
                               : ml_copysign (ml_make_inf (0), s);
        return (cplx) {re, im};
    }

    double c = ml_cos (a.imag);
    double s = ml_sin (a.imag);

    return (cplx) {mag * c, mag * s};
}
/* Long-double complex logarithm for internal consumers (power path).
 * Round-4: the public double log is 0 ULP, but at |log|~691 that is still
 * 5.7e-14 ABSOLUTE, and b*log amplifies it (3.7x here -> 885 ULP in pow).
 * Keeping (re,im) in LD preserves ~1e-19 and fixes the whole chain.
 * Same branch structure/edge behavior as ml_cplx_logarithm. */
static void ml_cplx_logarithm_ld (cplx a, long double *re, long double *im) {
    /* Zero base is handled by callers (a_zero early-return in pow paths);
     * guard here so a future reorder cannot turn log(0) into Inf*0=NaN
     * poisoning: re=-Inf, im=atan2(0,0)=0 by convention. */
    if (a.real == 0.0 && a.imag == 0.0) {
        *re = -__builtin_huge_vall ();
        *im = 0.0L;
        return;
    }
    long double xx = (long double) a.real, yy = (long double) a.imag;
    long double m = xx * xx + yy * yy;
    if (m == 0.0L) {
        *re = -__builtin_huge_vall ();
    } else if (m >= 0.5L && m <= 2.0L) {
        *re = 0.5L * __builtin_log1pl (m - 1.0L);
    } else {
        *re = 0.5L * __builtin_logl (m);
    }
    *im = (long double) ml_cplx_arg (a);
    /* arg in LD when it matters: ml_cplx_arg is 1 ULP double (~9e-17
     * absolute at pi/4); b*arg amplifies the same way. Use atan2l. */
    if (ml_isfinite (a.real) && ml_isfinite (a.imag) && !(a.real == 0.0 && a.imag == 0.0)) {
        *im = __builtin_atan2l (yy, xx);
    }
}
ML_API cplx ml_cplx_logarithm (cplx a) {
    /* Round-4: real part as 0.5*log(x^2+y^2) in long double (was
     * log(hypot): hypot's 0.35-ULP error at |a|~=1 amplified 88x by log
     * into log(1+0.1i).real (output 0.005 vs abs error 7.6e-17).
     * Even correctly-rounded hypot leaves 126 ULP here — the |a| must
     * carry extended precision. LD x^2+y^2 has 64-bit mantissa; m-1 is
     * Sterbenz-exact on [0.5,2] and log1pl is 1 ULP LD: total ~0.6 ULP
     * double. Outside [0.5,2] plain 0.5*logl (no amplification regime).
     * Edge behavior preserved (|a|=0 -> -inf; inf -> inf; NaN -> NaN). */
    if (ml_isnan (a.real) || ml_isnan (a.imag)) {
        double n = ml_make_nan ();
        return (cplx) {n, ml_cplx_arg (a)};
    }
#if defined(__STDC_VERSION__) && defined(__LDBL_MANT_DIG__) && (__LDBL_MANT_DIG__ >= 64)
    {
        long double xx = (long double) a.real, yy = (long double) a.imag;
        long double m = xx * xx + yy * yy;
        /* m == 0 (a == 0, since xx*xx+yy*yy underflows to 0 only for
         * subnormal-square inputs whose log is still finite... except true
         * a==0): route exact zero to -inf; any m>0 gives finite re here
         * (m <= 2*DBL_MAX^2 -> re <= 710.5) so a plain cast suffices. */
        if (m == 0.0L) {
            if (a.real == 0.0 && a.imag == 0.0)
                return (cplx) {-ml_make_inf (0), ml_cplx_arg (a)};
        }
        {
            long double re = (m >= 0.5L && m <= 2.0L)
                                 ? 0.5L * __builtin_log1pl (m - 1.0L)
                                 : 0.5L * __builtin_logl (m);
            return (cplx) {(double) re, ml_cplx_arg (a)};
        }
    }
#else
    return (cplx) {ml_log1p (ml_cplx_abs (a) - 1.0), ml_cplx_arg (a)};
#endif
}
ML_API cplx ml_cplx_sqrt (cplx a) {
    /* Hull-Beebe stable form (no cancellation):
     *   re(a)>=0: t=sqrt(2(r+x)); re=t/2; im=y/t
     *   re(a)<0:  t=sqrt(2(r-x)); re=|y|/t; im=sign(y)*t/2.
     * Round-4: the old im=sqrt((r-x)/2) cancels catastrophically when
     * |y|<<|x| (sqrt(1e300+1e200i) gave im=0 instead of 5e49: hypot
     * rounds 1+1e-200 to 1, so r-x underflows the true 1e100 gap to 0;
     * cmath gets 5e49 via the division form). The division form is exact
     * here: the *2 scaling is exact (power of two) and r+x is
     * well-conditioned (same-sign sum, no cancellation), sqrt 1 ULP. */
    if (ml_isnan (a.real) || ml_isnan (a.imag)) {
        double n = ml_make_nan ();
        return (cplx) {n, n};
    }
    if (ml_isinf (a.real) || ml_isinf (a.imag)) {
        /* Infinite magnitude dominates: check imag-Inf first so
         * (-Inf+Infi) gives (Inf,Inf), not (0,Inf). */
        if (ml_isinf (a.imag)) {
            if (ml_isinf (a.real)) {
                if (a.real > 0.0) {
                    return (cplx) {ml_make_inf (0), ml_copysign (ml_make_inf (0), a.imag)};
                }
                return (cplx) {ml_make_inf (0), ml_copysign (ml_make_inf (0), a.imag)};
            }
            return (cplx) {ml_make_inf (0), ml_copysign (ml_make_inf (0), a.real)};
        }
        if (ml_isinf (a.real)) {
            if (a.real > 0.0) {
                return (cplx) {ml_make_inf (0), ml_copysign (ml_make_inf (0), a.imag)};
            }
            return (cplx) {0.0, ml_copysign (ml_make_inf (0), a.imag)};
        }
        return (cplx) {ml_make_inf (0), ml_copysign (ml_make_inf (0), a.real)};
    }
    {
        double r = ml_cplx_abs (a);
        /* C99 Annex G: csqrt(-0,+0) = (+0,+0). Force +0 real for zero input. */
        if (r == 0.0) return (cplx) {0.0, ml_copysign (0.0, a.imag)};
        if (a.real >= 0.0) {
            double t = ml_sqrt (2.0 * (r + a.real));
            if (t == 0.0) return (cplx) {0.0, ml_copysign (0.0, a.imag)};
            return (cplx) {t * 0.5, a.imag / t};
        } else {
            double t = ml_sqrt (2.0 * (r - a.real));
            if (t == 0.0) return (cplx) {ml_fabs (a.imag), ml_copysign (0.0, a.imag)};
            {
                double re = ml_fabs (a.imag) / t;
                double im = ml_copysign (t * 0.5, a.imag);
                return (cplx) {re, im};
            }
        }
    }
}
/* MATHLIB_CLOSURE_P0_COMPLEX_POWER_TREE */
static int ml_cplx_is_nan (cplx z) {
    return ml_isnan (z.real) || ml_isnan (z.imag);
}
/* Long-double b*log(a) exponential for the power paths below.
 * Round-4: the old ml_cplx_exponential({b.re*l.re-b.im*l.im, ...})
 * rounded w in doubles, so direction/magnitude error scaled with |w|
 * (3.5e-13 at |w|~=2500). Computing w in LD and exponentiating in LD
 * (sinl/cosl resolve |wi|<=2000 to ~1 ULP) extends the clean range 100x.
 * Outside the LD window the double/Payne-Hanek fallback preserves the
 * existing (documented) behavior: garbage direction past ~1e19 needs
 * 65+ bit arithmetic no double implementation has (glibc agrees). */
static cplx ml_cplx_pow_exp_ld (long double wr, long double wi,
                                double fbr, double fbi) {
    /* LD direction error is ~|wi|/2200 ULP (vs |wi|/2 ULP double), so the
     * LD window extends to |wi|<=1e15; magnitude window covers LD exp
     * range with double-safe cast (overflow -> signed inf like the
     * fallback: m>0 so (double)(m*c) carries sign(c), matching the
     * copysign pattern, including c==0 -> signed zero). */
    if (wi > -1.0e15L && wi < 1.0e15L && wr > -11356.0L && wr < 11350.0L) {
        long double m = __builtin_expl (wr);
        long double c = __builtin_cosl (wi), s = __builtin_sinl (wi);
        double re = (double) (m * c), im = (double) (m * s);
        if (ml_isfinite (re) && ml_isfinite (im)) return (cplx) {re, im};
    }
    return ml_cplx_exponential ((cplx) {fbr, fbi});
}
/* Shared w = b*log(a) helper: LD pair for accuracy, double pair for the
 * outside-window fallback. Single source so the three pow paths cannot
 * drift (previously the LD product was replicated 3x). */
static void ml_cplx_pow_w (cplx b, long double lre, long double lim,
                           cplx log_a,
                           long double *wr, long double *wi,
                           double *fbr, double *fbi) {
    *wr = (long double) b.real * lre - (long double) b.imag * lim;
    *wi = (long double) b.real * lim + (long double) b.imag * lre;
    *fbr = b.real * log_a.real - b.imag * log_a.imag;
    *fbi = b.real * log_a.imag + b.imag * log_a.real;
}
static int ml_cplx_is_zero (cplx z) {
    return z.real == 0.0 && z.imag == 0.0;
}
static int ml_cplx_is_one (cplx z) {
    return z.real == 1.0 && z.imag == 0.0;
}
static int ml_cplx_is_pure_real_inf (cplx z) {
    return ml_isinf (z.real) && z.imag == 0.0;
}
static cplx ml_cplx_make_nan (void) {
    double n = ml_make_nan ();
    return (cplx) {n, n};
}
static cplx ml_cplx_make_zero (void) {
    return (cplx) {0.0, 0.0};
}
static cplx ml_cplx_make_one (void) {
    return (cplx) {1.0, 0.0};
}
static cplx ml_cplx_make_pos_inf (void) {
    return (cplx) {ml_make_inf (0), 0.0};
}
static cplx ml_cplx_make_neg_inf (void) {
    return (cplx) {-ml_make_inf (0), 0.0};
}
static cplx ml_cplx_make_neg_zero (void) {
    return (cplx) {ml_copysign (0.0, -1.0), 0.0};
}
ML_API cplx ml_cplx_power (cplx a, cplx b) {
    int a_nan = ml_cplx_is_nan (a);
    int b_nan = ml_cplx_is_nan (b);

    int a_zero = ml_cplx_is_zero (a);
    int a_one = ml_cplx_is_one (a);
    int b_zero = ml_cplx_is_zero (b);

    int a_inf = ml_isinf (a.real) || ml_isinf (a.imag);
    int b_inf = ml_isinf (b.real) || ml_isinf (b.imag);

    if (b_zero) {
        return ml_cplx_make_one ();
    }

    if (a_one) {
        return ml_cplx_make_one ();
    }

    if (a_nan) {
        return ml_cplx_make_nan ();
    }

    if (b_nan) {
        return ml_cplx_make_nan ();
    }

    if (a_zero) {
        if (b_inf) {
            if (ml_isinf (b.imag)) {
                return ml_cplx_make_nan ();
            }

            if (b.real > 0.0) {
                return ml_cplx_make_zero ();
            }

            if (b.real < 0.0) {
                return ml_cplx_make_pos_inf ();
            }

            return ml_cplx_make_nan ();
        }

        if (b.real > 0.0) {
            return ml_cplx_make_zero ();
        }

        if (b.real < 0.0) {
            return ml_cplx_make_pos_inf ();
        }

        return ml_cplx_make_nan ();
    }

    if (a_inf) {
        if (b_inf) {
            /* Pure-real infinite base and exponent: resolve by magnitude.
             * +Inf^+Inf=+Inf, +Inf^-Inf=+0; -Inf follows odd/even sign
             * when b is a finite odd integer, else magnitude rule. */
            if (ml_cplx_is_pure_real_inf (a) && b.imag == 0.0) {
                if (b.real > 0.0) {
                    if (a.real < 0.0 && ml_is_odd_integer_double (b.real)) {
                        return ml_cplx_make_neg_inf ();
                    }
                    return ml_cplx_make_pos_inf ();
                }
                if (b.real < 0.0) {
                    if (a.real < 0.0 && ml_is_odd_integer_double (b.real)) {
                        return ml_cplx_make_neg_zero ();
                    }
                    return ml_cplx_make_zero ();
                }
            }
            /* Non-pure-real infinities: fall through to exp(b*log a)
             * rather than unconditionally returning NaN. */
            if (!ml_cplx_is_pure_real_inf (a) || b.imag != 0.0) {
                long double lre, lim, wr, wi;
                double fbr, fbi;
                ml_cplx_logarithm_ld (a, &lre, &lim);
                cplx log_a = ml_cplx_logarithm (a);
                ml_cplx_pow_w (b, lre, lim, log_a, &wr, &wi, &fbr, &fbi);
                return ml_cplx_pow_exp_ld (wr, wi, fbr, fbi);
            }
            return ml_cplx_make_nan ();
        }

        if (!ml_cplx_is_pure_real_inf (a)) {
            return ml_cplx_make_nan ();
        }

        if (b.imag != 0.0) {
            return ml_cplx_make_nan ();
        }

        if (b.real > 0.0) {
            if (a.real < 0.0 && ml_is_odd_integer_double (b.real)) {
                return ml_cplx_make_neg_inf ();
            }
            return ml_cplx_make_pos_inf ();
        }

        if (b.real < 0.0) {
            if (a.real < 0.0 && ml_is_odd_integer_double (b.real)) {
                return ml_cplx_make_neg_zero ();
            }
            return ml_cplx_make_zero ();
        }

        return ml_cplx_make_one ();
    }

    if (b_inf) {
        /* Finite nonzero a with infinite b: resolve via exp(b*log a)
         * instead of unconditional NaN. |a|>1 -> Inf, |a|<1 -> 0,
         * |a|==1 -> indeterminate (NaN). Zero-guard the complex product:
         * IEEE Inf*0=NaN would poison the real-real case
         * (e.g. cpow(2,Inf+0i) imag must be 0, not NaN). */
        cplx log_a = ml_cplx_logarithm (a);
        if (!ml_isfinite (log_a.real) || !ml_isfinite (log_a.imag)) {
            return ml_cplx_make_nan ();
        }
        if (log_a.real == 0.0) {
            return ml_cplx_make_nan ();
        }
        {
            /* Round-4: w in LD from the LD logarithm (double log_a caps
             * absolute accuracy at 0.5 ULP of |log|, amplified by b).
             * lre/lim are always finite here (finite nonzero a); the
             * double re/im feed only the fallback path. Guard structure
             * unchanged (avoids Inf*0=NaN poisoning). */
            long double lre, lim;
            double re, im;
            long double wrl = 0.0L, wil = 0.0L;
            ml_cplx_logarithm_ld (a, &lre, &lim);
            if (b.imag == 0.0 && log_a.imag == 0.0) {
                re = b.real * log_a.real;
                im = 0.0;
                wrl = (long double) b.real * lre;
                wil = 0.0L;
            } else if (b.imag == 0.0) {
                re = b.real * log_a.real;
                im = b.real * log_a.imag;
                wrl = (long double) b.real * lre;
                wil = (long double) b.real * lim;
            } else if (log_a.imag == 0.0 && log_a.real != 0.0) {
                re = b.real * log_a.real;
                im = b.imag * log_a.real;
                wrl = (long double) b.real * lre;
                wil = (long double) b.imag * lre;
            } else {
                re = b.real * log_a.real - b.imag * log_a.imag;
                im = b.real * log_a.imag + b.imag * log_a.real;
                wrl = (long double) b.real * lre - (long double) b.imag * lim;
                wil = (long double) b.real * lim + (long double) b.imag * lre;
            }
            return ml_cplx_pow_exp_ld (wrl, wil, re, im);
        }
    }

    /* General path: w = b*log(a) entirely in LD (see above). */
    {
        long double lre, lim, wr, wi;
        double fbr, fbi;
        ml_cplx_logarithm_ld (a, &lre, &lim);
        cplx log_a = ml_cplx_logarithm (a);
        ml_cplx_pow_w (b, lre, lim, log_a, &wr, &wi, &fbr, &fbi);
        return ml_cplx_pow_exp_ld (wr, wi, fbr, fbi);
    }
}
ML_API cplx ml_cplx_nth_root (cplx a, int n, int k) {
    double nn = ml_make_nan ();
    if (ML_UNLIKELY (n <= 0)) return (cplx) {nn, nn};
    {
        int kk = ((k % n) + n) % n;
        double r = ml_cplx_abs (a);
        double th = ml_cplx_arg (a);
        if (!ml_isfinite (r) || ml_isnan (th)) {
            /* Infinite magnitude -> infinite root (angle moot).
             * Use signed Inf/zero to avoid Inf*0=NaN. */
            if (ml_isinf (r)) {
                double ph = 0.0;
                /* Preserve direction when defined. */
                if (ml_isfinite (th)) ph = (th + 2.0 * ML_PI * (double) kk) / (double) n;
                {
                    double c = ml_cos (ph), s = ml_sin (ph);
                    double re = (c == 0.0) ? ml_copysign (0.0, c) : ml_copysign (ml_make_inf (0), c);
                    double im = (s == 0.0) ? ml_copysign (0.0, s) : ml_copysign (ml_make_inf (0), s);
                    return (cplx) {re, im};
                }
            }
            return (cplx) {nn, nn};
        }
        if (r == 0.0) return (cplx) {0.0, 0.0};
        {
            double m = ml_exp (ml_log (r) / (double) n);
            double ph = (th + 2.0 * ML_PI * (double) kk) / (double) n;
            return (cplx) {m * ml_cos (ph), m * ml_sin (ph)};
        }
    }
}
