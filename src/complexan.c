#include "ml_compiler.h"
#include "ml_complexan.h"
#include "ml_trig.h"

/* ============================================================================
 * v12R2 complex analysis pack.
 * ========================================================================== */

#define CA_PI 3.14159265358979323846
static cplx ca_make (double re, double im) {
    cplx z;
    z.real = re;
    z.imag = im;
    return z;
}
static double ca_abs (cplx z) { return ml_sqrt (z.real * z.real + z.imag * z.imag); }
ML_API cplx ml_cplx_poly_eval (const cplx *c, int n, cplx z) {
    cplx r = ca_make (0.0, 0.0);
    int i;
    if (ML_UNLIKELY (!c || n < 0)) return r;
    for (i = 0; i <= n; i++)
        r = ml_cplx_add (ml_cplx_mul (r, z), c [i]);
    return r;
}
ML_API ml_status_t ml_cplx_poly_roots (const cplx *c, int n, cplx *roots) {
    int i, j, it;
    double lead;
    cplx norm [32];
    if (ML_UNLIKELY (!c || !roots || n <= 0 || n > 31)) return ML_ERR_INVALID_ARG;
    lead = ca_abs (c [0]);
    if (lead == 0.0) return ML_ERR_SINGULAR;
    for (i = 0; i <= n; i++)
        norm [i] = ml_cplx_div (c [i], ca_make (c [0].real, c [0].imag));
    for (i = 0; i < n; i++) {
        double ang = 2.0 * CA_PI * (double) i / (double) n;
        roots [i] = ca_make (0.4 * ml_cos (ang), 0.9 * 0.4 * ml_sin (ang));
        if (ca_abs (roots [i]) == 0.0) roots [i] = ca_make (0.1, 0.1);
    }
    for (it = 0; it < 200; it++) {
        double mx = 0.0;
        for (i = 0; i < n; i++) {
            cplx num = ml_cplx_poly_eval (norm, n, roots [i]);
            cplx den = ca_make (1.0, 0.0);
            cplx delta;
            for (j = 0; j < n; j++) {
                if (j == i) continue;
                den = ml_cplx_mul (den, ml_cplx_sub (roots [i], roots [j]));
            }
            if (ca_abs (den) < 1e-300) den = ca_make (1e-300, 1e-300);
            delta = ml_cplx_div (num, den);
            roots [i] = ml_cplx_sub (roots [i], delta);
            if (ca_abs (delta) > mx) mx = ca_abs (delta);
        }
        if (mx < 1e-14) break;
    }
    return ML_SUCCESS;
}
ML_API cplx ml_cplx_sin (cplx z) {
    return ca_make (ml_sin (z.real) * ml_cosh (z.imag),
                    ml_cos (z.real) * ml_sinh (z.imag));
}
ML_API cplx ml_cplx_cos (cplx z) {
    return ca_make (ml_cos (z.real) * ml_cosh (z.imag),
                    -ml_sin (z.real) * ml_sinh (z.imag));
}
ML_API cplx ml_cplx_tan (cplx z) {
    cplx s = ml_cplx_sin (z), c = ml_cplx_cos (z);
    return ml_cplx_div (s, c);
}
ML_API cplx ml_cplx_derivative (ml_cplx_func_t f, cplx z, void *ctx) {
    double h = 1e-6;
    cplx zp = ca_make (z.real + h, z.imag + h);
    cplx zm = ca_make (z.real - h, z.imag - h);
    cplx a = f (zp, ctx), b = f (zm, ctx);
    cplx d = ml_cplx_sub (a, b);
    cplx den = ca_make (2.0 * h, 2.0 * h);
    return ml_cplx_div (d, den);
}
ML_API cplx ml_cplx_circle_integral (ml_cplx_func_t f, cplx center,
                                     double radius, int n, void *ctx) {
    cplx sum = ca_make (0.0, 0.0);
    int k;
    if (ML_UNLIKELY (!f || n <= 0 || radius <= 0.0)) return sum;
    for (k = 0; k < n; k++) {
        double th = 2.0 * CA_PI * (double) k / (double) n;
        cplx z = ca_make (center.real + radius * ml_cos (th),
                          center.imag + radius * ml_sin (th));
        cplx dz = ca_make (-radius * ml_sin (th), radius * ml_cos (th));
        cplx term = ml_cplx_mul (f (z, ctx), dz);
        sum = ml_cplx_add (sum, term);
    }
    {
        double dth = 2.0 * CA_PI / (double) n;
        sum.real *= dth;
        sum.imag *= dth;
    }
    return sum;
}
ML_API cplx ml_cplx_residue (ml_cplx_func_t f, cplx pole, double radius, int n,
                             void *ctx) {
    cplx I = ml_cplx_circle_integral (f, pole, radius, n, ctx);
    /* divide by 2*pi*i:  (a+bi)/(2*pi*i) = (b - a*i)/(2*pi) */
    return ca_make (I.imag / (2.0 * CA_PI), -I.real / (2.0 * CA_PI));
}
ML_API void ml_cplx_solve_quadratic (cplx a, cplx b, cplx c, cplx *r1, cplx *r2) {
    cplx disc, sq, two_a, negb;
    if (ML_UNLIKELY (!r1 || !r2)) return;
    two_a = ca_make (2.0 * a.real, 2.0 * a.imag);
    negb = ca_make (-b.real, -b.imag);
    disc = ml_cplx_sub (ml_cplx_mul (b, b), ml_cplx_mul (ca_make (4.0, 0.0),
                                                         ml_cplx_mul (a, c)));
    sq = ml_cplx_sqrt (disc);
    *r1 = ml_cplx_div (ml_cplx_add (negb, sq), two_a);
    *r2 = ml_cplx_div (ml_cplx_sub (negb, sq), two_a);
}
ML_API int ml_cplx_escape (cplx c, int max_iter, double radius) {
    cplx z = ca_make (0.0, 0.0);
    int it;
    if (max_iter <= 0) max_iter = 100;
    if (radius <= 0.0) radius = 2.0;
    for (it = 0; it < max_iter; it++) {
        z = ml_cplx_add (ml_cplx_mul (z, z), c);
        if (ca_abs (z) > radius) return it;
    }
    return max_iter;
}
