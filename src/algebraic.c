#include "ml_compiler.h"
#include "ml_algebraic.h"
#include "ml_core.h"
#include "ml_trig.h"
#include "ml_topology.h"
#include <stddef.h>
#include <stdlib.h>

/* Resultants, discriminants, Sturm isolation, hook formula, Rodrigues. */

#define MLA_CMAX 1000000LL
#define MLA_DEGMAX 12
static int mla_valid_poly (const long long *p, int deg) {
    int i;
    if (!p || deg < 1 || deg > MLA_DEGMAX) return 0;
    if (p [deg] == 0) return 0;
    for (i = 0; i <= deg; i++) {
        if (p [i] > MLA_CMAX || p [i] < -MLA_CMAX) return 0;
    }
    return 1;
}
ML_API ml_status_t ml_resultant (const long long *p, int degp,
                                 const long long *q, int degq, long long *res) {
    /* Sylvester matrix (degp+degq) x (degp+degq), Bareiss det. */
    int n;
    long long *mat = NULL;
    int i;
    int j;
    ml_status_t st;
    if (ML_UNLIKELY (!p || !q || !res)) return ML_ERR_INVALID_ARG;
    if (!mla_valid_poly (p, degp) || !mla_valid_poly (q, degq)) return ML_ERR_INVALID_ARG;
    n = degp + degq;
    if ((size_t) n * (size_t) n > (size_t) -1 / sizeof (long long)) return ML_ERR_INVALID_ARG;
    mat = (long long *) malloc ((size_t) n * (size_t) n * sizeof (long long));
    if (!mat) return ML_ERR_WORKSPACE;
    for (i = 0; i < n * n; i++) mat [i] = 0;
    for (i = 0; i < degq; i++) {
        for (j = 0; j <= degp; j++) mat [(size_t) i * (size_t) n + (size_t) (i + j)] = p [j];
    }
    for (i = 0; i < degp; i++) {
        for (j = 0; j <= degq; j++) mat [(size_t) (degq + i) * (size_t) n + (size_t) (i + j)] = q [j];
    }
    st = ml_bareiss_det (mat, n, res);
    free (mat);
    return st;
}
ML_API ml_status_t ml_discriminant (const long long *p, int deg, long long *disc) {
    /* disc = (-1)^{n(n-1)/2} Res(p,p')/lead. */
    long long dp [12];
    long long res = 0;
    ml_status_t st;
    int i;
    if (ML_UNLIKELY (!p || !disc)) return ML_ERR_INVALID_ARG;
    /* Linear edge first: discriminant is 1 by convention (single root). */
    if (deg == 1 && p [0] <= MLA_CMAX && p [0] >= -MLA_CMAX &&
        p [1] != 0 && p [1] <= MLA_CMAX && p [1] >= -MLA_CMAX) {
        *disc = 1;
        return ML_SUCCESS;
    }
    if (!mla_valid_poly (p, deg)) return ML_ERR_INVALID_ARG;
    for (i = 0; i < deg; i++) dp [i] = p [i + 1] * (long long) (i + 1);
    st = ml_resultant (p, deg, dp, deg - 1, &res);
    if (st != ML_SUCCESS) return st;
    {
        long long lead = p [deg];
        long long num = res;
        int halves = (deg * (deg - 1) / 2) & 1;
        /* Divide by lead first (discriminant identity), then sign. */
        if (lead == 0) return ML_ERR_INTERNAL;
        if (num % lead != 0) return ML_ERR_INTERNAL;
        num /= lead;
        *disc = halves ? -num : num;
        return ML_SUCCESS;
    }
}
#define MLA_STURM_DEG 8
static int mla_sturm_sign (long double v, long double scale) {
    long double tol = 1e-9L * (1.0L + scale);
    if (v > tol) return 1;
    if (v < -tol) return -1;
    return 0;
}
static long double mla_poly_eval_ld (const long double *c, int deg, long double x) {
    long double r = 0.0L;
    int i;
    for (i = deg; i >= 0; i--) r = r * x + c [i];
    return r;
}
ML_API int ml_sturm_count (const long long *p, int deg, double a, double b) {
    /* Sturm chain by LD long division (degrees <= 8, small ints exact);
     * count = variations(a) - variations(b), roots in (a, b]. */
    long double seq [9][9];
    int lens [9];
    int nseq = 0;
    int i;
    int va = 0;
    int vb = 0;
    int prev = 0;
    if (!p || deg < 1 || deg > MLA_STURM_DEG) return -1;
    if (ml_isnan (a) || ml_isnan (b) || !ml_isfinite (a) || !ml_isfinite (b)) return -1;
    if (!(a < b)) return -1;
    for (i = 0; i <= deg; i++) {
        if (p [i] > MLA_CMAX || p [i] < -MLA_CMAX) return -1;
    }
    if (p [deg] == 0) return -1;
    for (i = 0; i <= 8; i++) {
        int j;
        for (j = 0; j <= 8; j++) seq [i][j] = 0.0L;
        lens [i] = 0;
    }
    for (i = 0; i <= deg; i++) seq [0][i] = (long double) p [i];
    lens [0] = deg;
    for (i = 0; i < deg; i++) seq [1][i] = (long double) p [i + 1] * (long double) (i + 1);
    lens [1] = deg - 1;
    nseq = 2;
    /* Reduce: seq[k] = -rem(seq[k-2], seq[k-1]). */
    while (nseq < 9) {
        long double *A = seq [nseq - 2];
        long double *B = seq [nseq - 1];
        int da = lens [nseq - 2];
        int db = lens [nseq - 1];
        long double R [9];
        int dr;
        int k;
        for (k = 0; k < 9; k++) R [k] = 0.0L;
        if (db < 0) break;
        if (db == 0 && mla_sturm_sign (B [0], 1.0L) == 0) break;
        if (db == 0) {
            /* Constant nonzero divisor: remainder is zero; chain ends. */
            break;
        }
        dr = da;
        {
            long double W [9];
            for (k = 0; k <= da; k++) W [k] = A [k];
            while (dr >= db && dr >= 0) {
                long double lead = B [db];
                long double coef;
                int e;
                if (mla_sturm_sign (lead, 1.0L) == 0) {
                    dr = -1;
                    break;
                }
                coef = W [dr] / lead;
                e = dr - db;
                for (k = 0; k <= db; k++) W [e + k] -= coef * B [k];
                while (dr >= 0 && mla_sturm_sign (W [dr], 1.0L) == 0) dr--;
            }
            for (k = 0; k <= dr && k < 9; k++) R [k] = -W [k];
        }
        while (dr >= 0 && mla_sturm_sign (R [dr], 1.0L) == 0) dr--;
        if (dr < 0) break;
        for (k = 0; k <= dr; k++) seq [nseq][k] = R [k];
        lens [nseq] = dr;
        nseq++;
        if (dr == 0) break;
    }
    /* Sign variations at a (skip zeros), minus variations at b. */
    prev = 0;
    for (i = 0; i < nseq; i++) {
        long double sc = 0.0L;
        int k;
        int s;
        for (k = 0; k <= lens [i]; k++) {
            long double v = seq [i][k];
            if (v < 0.0L) v = -v;
            if (v > sc) sc = v;
        }
        s = mla_sturm_sign (mla_poly_eval_ld (seq [i], lens [i], (long double) a), sc);
        if (s == 0) continue;
        if (prev != 0 && s != prev) va++;
        prev = s;
    }
    prev = 0;
    for (i = 0; i < nseq; i++) {
        long double sc = 0.0L;
        int k;
        int s;
        for (k = 0; k <= lens [i]; k++) {
            long double v = seq [i][k];
            if (v < 0.0L) v = -v;
            if (v > sc) sc = v;
        }
        s = mla_sturm_sign (mla_poly_eval_ld (seq [i], lens [i], (long double) b), sc);
        if (s == 0) continue;
        if (prev != 0 && s != prev) vb++;
        prev = s;
    }
    return va - vb;
}
ML_API int ml_descartes_bound (const long long *p, int deg) {
    /* Nonzero-coefficient sign changes (exact integer scan). */
    int changes = 0;
    int prev = 0;
    int i;
    if (!p || deg < 1 || deg > MLA_DEGMAX) return -1;
    for (i = 0; i <= deg; i++) {
        if (p [i] > MLA_CMAX || p [i] < -MLA_CMAX) return -1;
    }
    if (p [deg] == 0) return -1;
    for (i = 0; i <= deg; i++) {
        int s = (p [i] > 0) ? 1 : ((p [i] < 0) ? -1 : 0);
        if (s == 0) continue;
        if (prev != 0 && s != prev) changes++;
        prev = s;
    }
    return changes;
}
ML_API uint64_t ml_hook_dimension (const int *parts, int len) {
    /* Specht dim = n! / prod hooks, GCD-cancelled exact. n <= 20 keeps n!
     * in u64 (larger irreps overflow it genuinely). Ferrers rows required. */
    int n = 0;
    int i;
    int j;
    uint64_t num = 1;
    uint64_t den = 1;
    uint64_t g;
    if (!parts || len <= 0 || len > 20) return 0;
    for (i = 0; i < len; i++) {
        if (parts [i] <= 0 || parts [i] > 20) return 0;
        n += parts [i];
    }
    if (n > 20 || n <= 0) return 0;
    for (i = 1; i <= n; i++) {
        if (num > UINT64_MAX / (uint64_t) i) return UINT64_MAX;
        num *= (uint64_t) i;
    }
    for (i = 0; i < len; i++) {
        for (j = 0; j < parts [i]; j++) {
            /* Hook of cell (i,j): right + down + 1. Partition need not be
             * sorted for the multiset formula, but hooks assume Ferrers
             * rows; sort check: rows nonincreasing required. */
            int right;
            int down = 0;
            int k;
            uint64_t h;
            if (i > 0 && parts [i] > parts [i - 1]) return 0;
            right = parts [i] - j - 1;
            for (k = i + 1; k < len; k++) {
                if (parts [k] > j) down++;
            }
            h = (uint64_t) (right + down + 1);
            g = 0;
            {
                uint64_t a2 = num, b2 = h;
                while (b2 != 0) {
                    uint64_t t = a2 % b2;
                    a2 = b2;
                    b2 = t;
                }
                g = a2;
            }
            if (g > 1) {
                num /= g;
                h /= g;
            }
            g = 0;
            {
                uint64_t a2 = den, b2 = h;
                while (b2 != 0) {
                    uint64_t t = a2 % b2;
                    a2 = b2;
                    b2 = t;
                }
                g = a2;
            }
            if (g > 1) {
                den /= g;
                h /= g;
            }
            if (h != 1) {
                if (den > UINT64_MAX / h) return UINT64_MAX;
                den *= h;
            }
        }
    }
    if (den == 0 || num % den != 0) return UINT64_MAX;
    return num / den;
}
ML_API uint64_t ml_sn_class_size (const int *parts, int len) {
    /* n! / prod_k (m_k! k^{m_k}), m_k = multiplicity of part k. */
    int n = 0;
    int i;
    int k;
    uint64_t num = 1;
    if (!parts || len <= 0 || len > 60) return 0;
    for (i = 0; i < len; i++) {
        if (parts [i] <= 0 || parts [i] > 60) return 0;
        n += parts [i];
    }
    if (n > 20) return UINT64_MAX;
    for (i = 1; i <= n; i++) num *= (uint64_t) i;
    for (k = 1; k <= n; k++) {
        int m = 0;
        uint64_t kf = 1;
        uint64_t mf = 1;
        for (i = 0; i < len; i++) {
            if (parts [i] == k) m++;
        }
        if (m == 0) continue;
        for (i = 1; i <= m; i++) {
            if (mf > UINT64_MAX / (uint64_t) i) return UINT64_MAX;
            mf *= (uint64_t) i;
        }
        for (i = 0; i < m; i++) {
            if (kf > UINT64_MAX / (uint64_t) k) return UINT64_MAX;
            kf *= (uint64_t) k;
        }
        if (mf != 0 && kf > UINT64_MAX / mf) return UINT64_MAX;
        {
            uint64_t blk = mf * kf;
            if (blk == 0 || num % blk != 0) {
                /* Divide stepwise to keep exactness. */
                if (num % mf != 0) return UINT64_MAX;
                num /= mf;
                if (num % kf != 0) return UINT64_MAX;
                num /= kf;
            } else {
                num /= blk;
            }
        }
    }
    return num;
}
ML_API ml_status_t ml_rodrigues_exp (const double axis [3], double theta, double R [9]) {
    /* R = I + sin t K + (1-cos t) K^2; axis must be unit (tol 1e-9). */
    double nx;
    double ny;
    double nz;
    double s;
    double c;
    double norm;
    if (ML_UNLIKELY (!axis || !R)) return ML_ERR_INVALID_ARG;
    if (ml_isnan (theta) || !ml_isfinite (theta)) return ML_ERR_NAN_INPUT;
    nx = axis [0];
    ny = axis [1];
    nz = axis [2];
    if (ml_isnan (nx) || ml_isnan (ny) || ml_isnan (nz)) return ML_ERR_NAN_INPUT;
    if (!ml_isfinite (nx) || !ml_isfinite (ny) || !ml_isfinite (nz)) return ML_ERR_NAN_INPUT;
    norm = ml_sqrt (nx * nx + ny * ny + nz * nz);
    if (!(norm > 0.0) || !ml_isfinite (norm)) return ML_ERR_INVALID_ARG;
    if (ml_fabs (norm - 1.0) > 1e-9 * (1.0 + norm)) return ML_ERR_INVALID_ARG;
    nx /= norm;
    ny /= norm;
    nz /= norm;
    s = ml_sin (theta);
    c = ml_cos (theta);
    {
        double k [9] = {0.0, -nz, ny, nz, 0.0, -nx, -ny, nx, 0.0};
        double k2 [9];
        int i;
        int j;
        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                long double t2 = 0.0L;
                int l;
                for (l = 0; l < 3; l++) t2 += (long double) k [(size_t) i * 3 + (size_t) l] * (long double) k [(size_t) l * 3 + (size_t) j];
                k2 [(size_t) i * 3 + (size_t) j] = (double) t2;
            }
        }
        for (i = 0; i < 9; i++) {
            int r = i / 3;
            int cc = i % 3;
            double id = (r == cc) ? 1.0 : 0.0;
            R [i] = id + s * k [i] + (1.0 - c) * k2 [i];
            if (!ml_isfinite (R [i])) return ML_ERR_OVERFLOW;
        }
    }
    return ML_SUCCESS;
}
