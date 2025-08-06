#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <string.h>

// Define complex struct and PI constant for math clarity
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct {
    double real;
    double imag;
} Complex;


// --- START MERSENNE TWISTER (MT19937) --- Do Not Modify ---
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
// --- END MERSENNE TWISTER ---

// Global variables for benchmark data
static int num_samples;
static int num_transforms;
static Complex* input_signal;
static Complex* working_data;
static double final_result;


// In-place Radix-2 Decimation-In-Time FFT implementation
void fft_radix2(Complex* data, int n) {
    // 1. Bit-reversal permutation
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;

        if (i < j) {
            Complex temp = data[i];
            data[i] = data[j];
            data[j] = temp;
        }
    }

    // 2. Butterfly computations
    for (int len = 2; len <= n; len <<= 1) {
        double angle = -2.0 * M_PI / len;
        Complex wlen = {cos(angle), sin(angle)}; // Twiddle factor for this stage
        for (int i = 0; i < n; i += len) {
            Complex w = {1.0, 0.0}; // Twiddle factor for this butterfly
            for (int j = 0; j < len / 2; j++) {
                Complex u = data[i + j];
                Complex v = data[i + j + len / 2];

                // v * w
                Complex v_w;
                v_w.real = v.real * w.real - v.imag * w.imag;
                v_w.imag = v.real * w.imag + v.imag * w.real;

                // u + v*w & u - v*w
                data[i + j].real = u.real + v_w.real;
                data[i + j].imag = u.imag + v_w.imag;
                data[i + j + len / 2].real = u.real - v_w.real;
                data[i + j + len / 2].imag = u.imag - v_w.imag;

                // Update twiddle factor: w = w * wlen
                double w_temp_real = w.real * wlen.real - w.imag * wlen.imag;
                w.imag = w.real * wlen.imag + w.imag * wlen.real;
                w.real = w_temp_real;
            }
        }
    }
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_samples> <num_transforms_to_run> <seed>\n", argv[0]);
        exit(1);
    }

    num_samples = atoi(argv[1]);
    num_transforms = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    // Radix-2 FFT requires num_samples to be a power of 2
    if ((num_samples <= 0) || (num_samples & (num_samples - 1)) != 0) {
        fprintf(stderr, "FATAL: num_samples must be a positive power of 2.\n");
        exit(1);
    }

    mt_seed(seed);

    input_signal = (Complex*)malloc(num_samples * sizeof(Complex));
    working_data = (Complex*)malloc(num_samples * sizeof(Complex));

    if (!input_signal || !working_data) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Generate a single random input signal
    for (int i = 0; i < num_samples; i++) {
        input_signal[i].real = ((double)mt_rand() / (double)UINT32_MAX) * 2.0 - 1.0;
        input_signal[i].imag = ((double)mt_rand() / (double)UINT32_MAX) * 2.0 - 1.0;
    }

    final_result = 0.0;
}

void run_computation() {
    for (int i = 0; i < num_transforms; i++) {
        // Copy input to working buffer for each in-place transform
        memcpy(working_data, input_signal, num_samples * sizeof(Complex));

        // Perform the FFT
        fft_radix2(working_data, num_samples);

        // Accumulate a result to prevent dead code elimination.
        // We sum the magnitude of the first non-DC frequency bin (k=1).
        if (num_samples > 1) {
            Complex c = working_data[1];
            final_result += sqrt(c.real * c.real + c.imag * c.imag);
        }
    }
}

void cleanup() {
    free(input_signal);
    free(working_data);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final accumulated result to stdout
    printf("%.6f\n", final_result);

    // Print time_taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
