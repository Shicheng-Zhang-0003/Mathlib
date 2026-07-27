#include "test_harness.h"
#include "fft.h"
#include "ml_complex.h"

int main() {
    ml_test_ctx_t ctx;
    ml_test_init(&ctx, "DSP & FFT");
    printf("=== DSP & FFT Tests ===\n");

    cplx sig[4] = {{1,0}, {0,0}, {0,0}, {0,0}};
    ml_fft_execute(sig, 4);

    ASSERT_NEAR(&ctx, sig[0].real, 1.0, 1e-9, "FFT DC bin");
    ASSERT_NEAR(&ctx, sig[1].real, 1.0, 1e-9, "FFT Nyquist bin");

    return ml_test_summary(&ctx);
}
