#include <stdio.h>
#include <math.h>
#include <time.h>

// Capture raw libm functions BEFORE the Trojan Horse overrides them
static double (*raw_libm_sin)(double) = sin;
static double (*raw_libm_sqrt)(double) = sqrt;

#include "compat.h"
#include "simd_batch.h"

static inline unsigned long long rdtsc() {
    unsigned int lo, hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((unsigned long long)hi << 32) | lo;
}

volatile double sink = 0.0;

int main() {
    printf("=========================================================\n");
    printf("   MATHLIB v9: ANTI-LIBM SILICON BENCHMARK (rdtsc)\n");
    printf("=========================================================\n\n");

    #if defined(MATHLIB_PROFILE_GRAPHICS)
    printf("Active Profile: GRAPHICS (Minimax + Fast RSqrt)\n");
    #elif defined(MATHLIB_PROFILE_EMBEDDED)
    printf("Active Profile: EMBEDDED (CORDIC)\n");
    #else
    printf("Active Profile: SCIENTIFIC (IEEE-754 Strict)\n");
    #endif
    printf("---------------------------------------------------------\n\n");

    int iterations = 10000000;
    unsigned long long start, end;
    double cycles_per_op;

    start = rdtsc();
    for(int i = 0; i < iterations; i++) sink += raw_libm_sin(i * 0.0001);
    end = rdtsc();
    cycles_per_op = (double)(end - start) / iterations;
    printf("GNU libm sin()        : %6.2f cycles / op\n", cycles_per_op);

    start = rdtsc();
    for(int i = 0; i < iterations; i++) sink += ml_sin(i * 0.0001);
    end = rdtsc();
    cycles_per_op = (double)(end - start) / iterations;
    printf("mathlib ml_sin()      : %6.2f cycles / op\n", cycles_per_op);

    start = rdtsc();
    for(int i = 0; i < iterations; i++) sink += ml_sqrt(i + 1.0);
    end = rdtsc();
    cycles_per_op = (double)(end - start) / iterations;
    printf("mathlib ml_sqrt()     : %6.2f cycles / op\n", cycles_per_op);

    start = rdtsc();
    for(int i = 0; i < iterations; i++) sink += raw_libm_sqrt(i + 1.0);
    end = rdtsc();
    cycles_per_op = (double)(end - start) / iterations;
    printf("GNU libm sqrt()       : %6.2f cycles / op\n", cycles_per_op);

    double simd_in[4096] __attribute__((aligned(32)));
    double simd_out[4096] __attribute__((aligned(32)));
    for(int i=0; i<4096; i++) simd_in[i] = i + 1.0;

    start = rdtsc();
    for(int i = 0; i < 10000; i++) {
        for(int j=0; j<4096; j+=4) ml_simd_batch_poly(&simd_in[j], &simd_out[j]);
    }
    end = rdtsc();
    printf("mathlib SIMD Batch    : %6.2f cycles / 4 ops\n", (double)(end - start) / (10000.0 * 1024.0));

    printf("\n=========================================================\n");
    return 0;
}
