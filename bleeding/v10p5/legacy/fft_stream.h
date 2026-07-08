#ifndef MATHLIB_FFT_STREAM_H
#define MATHLIB_FFT_STREAM_H

#include "fft.h"
#include <stdlib.h>

typedef struct {
    cplx* buffer;
    int size;
    int write_idx;
} ml_fft_stream;

static inline ml_fft_stream ml_fft_stream_create(int size) {
    ml_fft_stream s;
    s.size = size;
    s.buffer = (cplx*)calloc(size, sizeof(cplx));
    s.write_idx = 0;
    return s;
}

static inline void ml_fft_stream_push(ml_fft_stream* s, double real_val) {
    s->buffer[s->write_idx].real = real_val;
    s->buffer[s->write_idx].imag = 0.0;
    s->write_idx = (s->write_idx + 1) % s->size;
}

static inline void ml_fft_stream_process(ml_fft_stream* s, cplx* out_spectrum) {
    for(int i=0; i<s->size; i++) {
        out_spectrum[i] = s->buffer[(s->write_idx + i) % s->size];
    }
    fft_execute(out_spectrum, s->size);
}

static inline void ml_fft_stream_free(ml_fft_stream* s) {
    free(s->buffer);
}
#endif
