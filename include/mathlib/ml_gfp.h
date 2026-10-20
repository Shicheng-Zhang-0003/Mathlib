#ifndef MATHLIB_ML_GFP_H
#define MATHLIB_ML_GFP_H

#include "ml_compiler.h"
#include "ml_types.h"
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
/* Prime field F_p for prime p <= UINT32_MAX. Elements use canonical residues
 * [0,p). Extension fields F_(p^n), polynomials over F_p, and cryptographic
 * constant-time guarantees are outside this API's current scope. */
typedef struct {
    uint32_t p;
    uint32_t initialized;
} ml_gfp_ctx_t;
ML_API void ml_gfp_ctx_init (ml_gfp_ctx_t *ctx);
ML_API ml_status_t ml_gfp_ctx_set_prime (ml_gfp_ctx_t *ctx, uint32_t p);
ML_API int ml_gfp_ctx_is_valid (const ml_gfp_ctx_t *ctx);
ML_API ml_status_t ml_gfp_from_i64 (const ml_gfp_ctx_t *ctx, int64_t value,
                                    uint32_t *out);
ML_API ml_status_t ml_gfp_add (const ml_gfp_ctx_t *ctx, uint32_t a, uint32_t b,
                               uint32_t *out);
ML_API ml_status_t ml_gfp_sub (const ml_gfp_ctx_t *ctx, uint32_t a, uint32_t b,
                               uint32_t *out);
ML_API ml_status_t ml_gfp_mul (const ml_gfp_ctx_t *ctx, uint32_t a, uint32_t b,
                               uint32_t *out);
ML_API ml_status_t ml_gfp_pow (const ml_gfp_ctx_t *ctx, uint32_t a,
                               uint64_t exponent, uint32_t *out);
ML_API ml_status_t ml_gfp_inv (const ml_gfp_ctx_t *ctx, uint32_t a, uint32_t *out);
ML_API ml_status_t ml_gfp_div (const ml_gfp_ctx_t *ctx, uint32_t a, uint32_t b,
                               uint32_t *out);
/* Returns -1, 0, or +1 (Euler criterion). */
ML_API ml_status_t ml_gfp_legendre (const ml_gfp_ctx_t *ctx, uint32_t a,
                                    int *symbol);
/* Principal root is the smaller of r and p-r. Nonsquares return SINGULAR. */
ML_API ml_status_t ml_gfp_sqrt (const ml_gfp_ctx_t *ctx, uint32_t a,
                                uint32_t *root);

#ifdef __cplusplus
}
#endif
#endif /* MATHLIB_ML_GFP_H */
