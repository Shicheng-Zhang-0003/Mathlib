#include "ml_compiler.h"
#include "ml_manifold.h"
ML_API ml_status_t ml_proj_sphere(const double *x, double *y, int n) {
    if (!x || !y || n <= 0 || n > 256) return ML_ERR_INVALID_ARG;
    long double s = 0;
    for (int i = 0; i < n; i++) { if (!ml_isfinite(x[i])) return ML_ERR_NAN_INPUT; s += (long double)x[i]*x[i]; }
    long double nm = __builtin_sqrtl(s);
    if (!(nm > 0)) return ML_ERR_SINGULAR;
    for (int i = 0; i < n; i++) y[i] = (double)((long double)x[i]/nm);
    return ML_SUCCESS;
}
ML_API ml_status_t ml_proj_stiefel(const double *A, double *Q, int m, int n) {
    if (!A || !Q || m <= 0 || n <= 0 || m < n || m > 32 || n > 32) return ML_ERR_INVALID_ARG;
    /* Thin-QR projection via MGS with one re-orthogonalization pass.
     * This orthonormalizes the columns; it is NOT the polar-factor
     * nearest-point projection, and rank-deficient A fails SINGULAR. */
    double V[32][32];
    for (int i = 0; i < m; i++) for (int j = 0; j < n; j++) { if (!ml_isfinite(A[i*n+j])) return ML_ERR_NAN_INPUT; V[i][j] = A[i*n+j]; }
    for (int k = 0; k < n; k++) {
        long double s = 0; for (int i = 0; i < m; i++) s += (long double)V[i][k]*V[i][k];
        long double nm = __builtin_sqrtl(s);
        if (!(nm > 0) || !ml_isfinite((double)nm)) return ML_ERR_SINGULAR;
        for (int i = 0; i < m; i++) V[i][k] = (double)((long double)V[i][k]/nm);
        for (int j = k + 1; j < n; j++) {
            long double d = 0; for (int i = 0; i < m; i++) d += (long double)V[i][k]*V[i][j];
            if (!ml_isfinite((double)d)) return ML_ERR_SINGULAR;
            for (int i = 0; i < m; i++) V[i][j] -= (double)(d * (long double)V[i][k]);
        }
    }
    /* One MGS re-pass for ill-conditioned A: re-project each column
     * against all previous ones and renormalize. */
    for (int k = 0; k < n; k++) {
        for (int i = 0; i < k; i++) {
            long double d = 0; for (int r = 0; r < m; r++) d += (long double)V[r][i]*V[r][k];
            if (!ml_isfinite((double)d)) return ML_ERR_SINGULAR;
            for (int r = 0; r < m; r++) V[r][k] -= (double)(d * (long double)V[r][i]);
        }
        {
            long double s = 0; for (int r = 0; r < m; r++) s += (long double)V[r][k]*V[r][k];
            long double nm = __builtin_sqrtl(s);
            if (!(nm > 0) || !ml_isfinite((double)nm)) return ML_ERR_SINGULAR;
            for (int r = 0; r < m; r++) V[r][k] = (double)((long double)V[r][k]/nm);
        }
    }
    for (int i = 0; i < m; i++) for (int j = 0; j < n; j++) Q[i*n+j] = V[i][j];
    return ML_SUCCESS;
}
ML_API ml_status_t ml_exp_sphere(const double *x, const double *v, double *y, int n) {
    if (!x || !v || !y || n <= 0 || n > 256) return ML_ERR_INVALID_ARG;
    long double nv = 0, dot = 0, nx = 0;
    for (int i = 0; i < n; i++) {
        if (!ml_isfinite(x[i]) || !ml_isfinite(v[i])) return ML_ERR_NAN_INPUT;
        nv += (long double)v[i]*v[i]; dot += (long double)x[i]*v[i];
        nx += (long double)x[i]*x[i];
    }
    /* DESPOT-AUDIT: enforce unit basepoint |x|=1 (was unchecked, silently
     * producing off-manifold y). Tolerance 1e-12 absolute on ||x|-1| plus
     * scale-aware tangency |<x,v>| <= 1e-12*(1+|x||v|). */
    {
        long double nmx = __builtin_sqrtl(nx);
        if (!ml_isfinite((double)nmx) || __builtin_fabsl(nmx - 1.0L) > 1e-12L)
            return ML_ERR_INVALID_ARG;
    }
    {
        long double sv = __builtin_sqrtl(nv);
        long double nmx = __builtin_sqrtl(nx);
        long double tol = 1e-12L * (1.0L + nmx * sv);
        if (__builtin_fabsl(dot) > (double)tol) return ML_ERR_INVALID_ARG;
    }
    long double th = __builtin_sqrtl(nv);
    if (th == 0) { for (int i = 0; i < n; i++) y[i] = x[i]; return ML_SUCCESS; }
    long double c = __builtin_cosl(th), s = __builtin_sinl(th) / th;
    for (int i = 0; i < n; i++) y[i] = (double)((long double)x[i]*c + (long double)v[i]*s);
    return ML_SUCCESS;
}
ML_API double ml_sphere_dist(const double *x, const double *y, int n) {
    if (!x || !y || n <= 0 || n > 4096) return ml_make_nan();
    /* DESPOT-AUDIT: stable small-angle form d=2*asin(|x-y|/2) instead of
     * acos(<x,y>) (cancellation for d~0). Falls back to acos-clamp only if
     * the chord path is non-finite. Inputs must be unit (tolerance 1e-9 on
     * |norm-1|); non-unit returns NaN rather than a misleading angle. */
    long double chord2 = 0, d = 0, nx = 0, ny = 0;
    for (int i = 0; i < n; i++) {
        if (!ml_isfinite(x[i]) || !ml_isfinite(y[i])) return ml_make_nan();
        long double df = (long double)x[i] - (long double)y[i];
        chord2 += df * df;
        d += (long double)x[i]*(long double)y[i];
        nx += (long double)x[i]*(long double)x[i];
        ny += (long double)y[i]*(long double)y[i];
    }
    {
        long double nmx = __builtin_sqrtl(nx), nmy = __builtin_sqrtl(ny);
        if (!ml_isfinite((double)nmx) || !ml_isfinite((double)nmy)) return ml_make_nan();
        if (__builtin_fabsl(nmx - 1.0L) > 1e-9L || __builtin_fabsl(nmy - 1.0L) > 1e-9L) return ml_make_nan();
    }
    {
        long double half = __builtin_sqrtl(chord2) * 0.5L;
        if (half < 0) half = 0;
        if (half > 1.0L) half = 1.0L;
        /* 2*asin is accurate for small d; for antipodal d~pi both forms
         * agree to ~1 ULP. */
        long double r = 2.0L * __builtin_asinl(half);
        double o = (double)r;
        if (ml_isfinite(o)) return o;
    }
    if (d > 1.0L) { d = 1.0L; }
    if (d < -1.0L) { d = -1.0L; }
    double r = (double)__builtin_acosl(d);
    return ml_isfinite(r) ? r : ml_make_nan();
}
