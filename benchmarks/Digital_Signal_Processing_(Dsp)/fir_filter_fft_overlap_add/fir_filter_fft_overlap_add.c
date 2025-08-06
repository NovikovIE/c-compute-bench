#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <complex.h>

// M_PI is non-standard, so define it if not available
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Mersenne Twister (MT19937) Generator
#define MT_N 624
#define MT_M 397
#define MT_MATRIX_A 0x9908b0dfUL
#define MT_UPPER_MASK 0x80000000UL
#define MT_LOWER_MASK 0x7fffffffUL

static uint32_t mt[MT_N];
static int mt_index = MT_N + 1;

void mt_seed(uint32_t seed) {
    mt[0] = seed;
    for (mt_index = 1; mt_index < MT_N; mt_index++) {
        mt[mt_index] = (1812433253UL * (mt[mt_index - 1] ^ (mt[mt_index - 1] >> 30)) + mt_index);
    }
}

uint32_t mt_rand(void) {
    uint32_t y;
    static const uint32_t mag01[2] = {0x0UL, MT_MATRIX_A};
    if (mt_index >= MT_N) {
        if (mt_index > MT_N) {
             fprintf(stderr, "FATAL: Mersenne Twister not seeded.\n");
             exit(1);
        }
        for (int i = 0; i < MT_N - MT_M; i++) {
            y = (mt[i] & MT_UPPER_MASK) | (mt[i + 1] & MT_LOWER_MASK);
            mt[i] = mt[i + MT_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (int i = MT_N - MT_M; i < MT_N - 1; i++) {
            y = (mt[i] & MT_UPPER_MASK) | (mt[i + 1] & MT_LOWER_MASK);
            mt[i] = mt[i + (MT_M - MT_N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[MT_N - 1] & MT_UPPER_MASK) | (mt[0] & MT_LOWER_MASK);
        mt[MT_N - 1] = mt[MT_M - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];
        mt_index = 0;
    }
    y = mt[mt_index++];
    y ^= (y >> 11); y ^= (y << 7) & 0x9d2c5680UL; y ^= (y << 15) & 0xefc60000UL; y ^= (y >> 18);
    return y;
}

// Type alias for complex numbers
typedef double complex CPLX;

// Global struct to hold benchmark data
struct BenchmarkData {
    int signal_length;
    int num_filter_taps;
    int fft_size;

    double* signal;
    double* filter_kernel;
    double* output_signal;

    // Pre-computed FFT of the filter
    CPLX* filter_fft;

    // Final result for validation
    double final_result;
} g_data;

// --- FFT Helper Functions ---

// In-place iterative Cooley-Tukey FFT
void fft_iterative(CPLX buf[], int n) {
    int log2n = log2(n);

    // Bit-reversal permutation
    for (int i = 0; i < n; i++) {
        int rev = 0;
        for (int j = 0; j < log2n; j++) {
            if ((i >> j) & 1) rev |= 1 << (log2n - 1 - j);
        }
        if (i < rev) {
            CPLX temp = buf[i];
            buf[i] = buf[rev];
            buf[rev] = temp;
        }
    }

    // Cooley-Tukey iterations
    for (int s = 1; s <= log2n; s++) {
        int m = 1 << s;
        int m2 = m >> 1;
        CPLX w_m = cexp(-I * M_PI / m2);
        for (int j = 0; j < n; j += m) {
            CPLX w = 1.0;
            for (int k = 0; k < m2; k++) {
                CPLX t = w * buf[j + k + m2];
                CPLX u = buf[j + k];
                buf[j + k] = u + t;
                buf[j + k + m2] = u - t;
                w *= w_m;
            }
        }
    }
}

// In-place IFFT using the FFT function
void ifft_iterative(CPLX buf[], int n) {
    for (int i = 0; i < n; i++) {
        buf[i] = conj(buf[i]);
    }
    fft_iterative(buf, n);
    for (int i = 0; i < n; i++) {
        buf[i] = conj(buf[i]) / n;
    }
}


void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s signal_length num_filter_taps fft_size seed\n", argv[0]);
        exit(1);
    }
    g_data.signal_length = atoi(argv[1]);
    g_data.num_filter_taps = atoi(argv[2]);
    g_data.fft_size = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    mt_seed(seed);

    if (g_data.fft_size <= 0 || (g_data.fft_size & (g_data.fft_size - 1)) != 0) {
        fprintf(stderr, "FATAL: fft_size must be a positive power of 2.\n");
        exit(1);
    }
    if (g_data.num_filter_taps >= g_data.fft_size) {
        fprintf(stderr, "FATAL: num_filter_taps must be smaller than fft_size.\n");
        exit(1);
    }

    int output_len = g_data.signal_length + g_data.num_filter_taps - 1;
    g_data.signal = (double*)malloc(g_data.signal_length * sizeof(double));
    g_data.filter_kernel = (double*)malloc(g_data.num_filter_taps * sizeof(double));
    g_data.output_signal = (double*)malloc(output_len * sizeof(double));
    g_data.filter_fft = (CPLX*)malloc(g_data.fft_size * sizeof(CPLX));

    if (!g_data.signal || !g_data.filter_kernel || !g_data.output_signal || !g_data.filter_fft) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    for (int i = 0; i < g_data.signal_length; i++) {
        g_data.signal[i] = ((double)mt_rand() / (double)UINT32_MAX) * 2.0 - 1.0;
    }
    for (int i = 0; i < g_data.num_filter_taps; i++) {
        g_data.filter_kernel[i] = ((double)mt_rand() / (double)UINT32_MAX);
    }
    for (int i = 0; i < output_len; i++) {
        g_data.output_signal[i] = 0.0;
    }

    // Pre-compute FFT of the zero-padded filter
    for (int i = 0; i < g_data.num_filter_taps; i++) {
        g_data.filter_fft[i] = g_data.filter_kernel[i];
    }
    for (int i = g_data.num_filter_taps; i < g_data.fft_size; i++) {
        g_data.filter_fft[i] = 0.0;
    }
    fft_iterative(g_data.filter_fft, g_data.fft_size);

    g_data.final_result = 0.0;
}

void run_computation() {
    int block_size = g_data.fft_size - g_data.num_filter_taps + 1;
    int output_len = g_data.signal_length + g_data.num_filter_taps - 1;

    CPLX* work_buffer = (CPLX*)malloc(g_data.fft_size * sizeof(CPLX));
    if (!work_buffer) {
        fprintf(stderr, "FATAL: Memory allocation failed in computation.\n");
        exit(1);
    }

    for (int i = 0; i < g_data.signal_length; i += block_size) {
        int current_block_size = g_data.signal_length - i;
        if (current_block_size > block_size) {
            current_block_size = block_size;
        }

        for (int k = 0; k < current_block_size; k++) {
            work_buffer[k] = g_data.signal[i + k];
        }
        for (int k = current_block_size; k < g_data.fft_size; k++) {
            work_buffer[k] = 0.0;
        }

        fft_iterative(work_buffer, g_data.fft_size);

        for (int k = 0; k < g_data.fft_size; k++) {
            work_buffer[k] *= g_data.filter_fft[k];
        }

        ifft_iterative(work_buffer, g_data.fft_size);

        for (int k = 0; k < g_data.fft_size; k++) {
            if (i + k < output_len) {
                g_data.output_signal[i + k] += creal(work_buffer[k]);
            }
        }
    }

    free(work_buffer);

    double sum = 0.0;
    for (int i = 0; i < output_len; i++) {
        sum += g_data.output_signal[i];
    }
    g_data.final_result = sum;
}

void cleanup() {
    free(g_data.signal);
    free(g_data.filter_kernel);
    free(g_data.output_signal);
    free(g_data.filter_fft);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%f\n", g_data.final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
