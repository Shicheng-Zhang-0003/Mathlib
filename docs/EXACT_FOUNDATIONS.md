# Exact and certified foundations

This note describes the first shared exact-arithmetic layer added to MathLib.
It is deliberately smaller than a computer algebra system or proof assistant.

## Ownership and errors

`ml_mpz_t`, `ml_mpq_t`, and `ml_polyq_t` own dynamically allocated storage.
Call the matching `*_init` before use and `*_clear` exactly once when done.
Do not copy these structs with assignment or `memcpy`; use their copy APIs.
All output objects must be initialized. Except for comparison/query functions,
operations build a temporary result and replace the output only on success.

Allocation failure returns `ML_ERR_WORKSPACE`; unrepresentable allocation
sizes return `ML_ERR_OVERFLOW`. Dynamic means limited by address space and
available memory, not literally unbounded. Callers should cap input sizes for
untrusted data.

## Dynamic integers and rationals

`ml_mpz_t` stores signed base-2^32 limbs. Addition, subtraction, multiplication,
division/remainder, shifts, GCD, and string conversion are exact. Multiplication
uses the schoolbook algorithm; division uses bitwise long division, so very
large inputs can be slow. Signed division truncates toward zero and the
remainder has the dividend's sign.

`ml_mpq_t` stores a reduced numerator and a positive denominator. Arithmetic,
comparison, and normalization are exact. `ml_mpq_to_interval` converts an exact
rational to either one exactly representable binary64 point or the two
adjacent binary64 values bracketing it. It locates the enclosure by exact
integer cross-multiplication, not by casting through `long double`.
`ml_mpq_set_decimal` parses finite decimal literals (including a base-10
exponent) as exact rationals; it accepts up to 10,000 significand digits and an
exponent magnitude of 10,000. Recurring decimals need an explicit rational or
an approximation policy.

The legacy `ml_bigint_t` remains a separate fixed-size 4096-bit type. It is
retained for its zero-allocation ABI and existing modular/probable-prime APIs.

## Outward-rounded intervals

`ml_interval_t` supports finite-endpoint addition, subtraction,
multiplication, division, negation, and hull. Arithmetic endpoints are widened
by one representable value after each nearest-rounded binary64 operation. The
arithmetic functions require `ml_interval_is_supported()` to report IEEE
binary64 and `FE_TONEAREST`; intervals that contain zero cannot be divisors.
If a finite enclosure cannot be represented, the operation returns an error
and leaves the output unchanged.

This contract does not cover transcendental functions. C/libm does not promise
a universal error bound for `exp`, `log`, trigonometric functions, or their
compositions, so those functions must not be added to this API until a
separate enclosure method and validation contract exist.

## Exact polynomials over Q

`ml_polyq_t` represents dense univariate polynomials in ascending coefficient
order. It provides exact arithmetic, differentiation, Horner evaluation,
Euclidean division, and monic GCD. The zero polynomial has no active
coefficients and degree -1.

`ml_polyq_verify_division` checks a proposed witness by recomputing the exact
identity `A = B*Q + R` and checking `degree(R) < degree(B)`. For nonzero B,
these conditions uniquely characterize quotient and remainder over Q[x]. This
is a narrow mathematical certificate verifier; it does not verify arbitrary
proofs or the correctness of every algorithm in the library.

The representation is dense and univariate. Sparse/multivariate polynomial
algorithms, expression DAGs, rewrite systems, factorization over Q, Gröbner
bases, and algebraic-number representations remain out of scope for this first
increment.

## Prime fields

`ml_gfp_ctx_t` represents F_p for a validated prime `p <= UINT32_MAX`.
Elements are canonical residues in `[0,p)`. Context initialization verifies
primality by trial division; arithmetic uses exact 64-bit intermediates.
Square roots use the p=2, p=3 mod 4, and Tonelli–Shanks cases, and return the
smaller representative of the two roots.

Only prime fields are implemented. Extension fields F_(p^n), polynomial rings
over fields, factorization, and constant-time/cryptographic guarantees are
not provided.

## Validation

`tests/test_edge_foundations.c` checks integers beyond 4096 bits, randomized
native-range arithmetic, large multiplication/division identities, rational
normalization, exact rational enclosures, and interval failure atomicity.
`tests/test_edge_polyq.c` checks Q[x] arithmetic and division certificates.
`tests/test_edge_gfp.c` checks prime validation and finite-field operations,
including Tonelli–Shanks and the largest 32-bit prime. These are regression
checks, not a machine-checked proof of the C implementation.
