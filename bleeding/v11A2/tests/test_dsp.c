#include "test_harness.h"
#include "fft.h"
#include "ml_complex.h"

int main() {
    printf("=== DSP & FFT Tests ===\n");
    cplx sig[4] = {{1,0}, {0,0}, {0,0}, {0,0}};
    fft_execute(sig, 4);
    ASSERT_NEAR(sig[0].real, 1.0, 1e-9, "FFT DC bin");
    ASSERT_NEAR(sig[1].real, 1.0, 1e-9, "FFT Nyquist bin");
    return test_harness_summary("DSP");
}