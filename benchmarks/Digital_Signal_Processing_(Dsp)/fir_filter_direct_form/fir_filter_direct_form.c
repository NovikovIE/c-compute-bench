#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- DO NOT MODIFY --- MERSENNE TWISTER GENERATOR ---
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
             fprintf(stderr, "FATAL: Mersenne Twister not seeded.");
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
// --- END OF MERSENNE TWISTER ---

// Benchmark parameters
static int signal_length;
static int num_filter_taps;

// Data arrays
static double* input_signal;
static double* filter_taps;
static double* output_signal;

// Result for verification/DCE prevention
static double accumulated_result;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <signal_length> <num_filter_taps> <seed>\n", argv[0]);
        exit(1);
    }

    signal_length = atoi(argv[1]);
    num_filter_taps = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (signal_length <= 0 || num_filter_taps <= 0) {
        fprintf(stderr, "ERROR: signal_length and num_filter_taps must be positive integers.\n");
        exit(1);
    }

    mt_seed(seed);

    input_signal = (double*)malloc(signal_length * sizeof(double));
    filter_taps = (double*)malloc(num_filter_taps * sizeof(double));
    output_signal = (double*)malloc(signal_length * sizeof(double));

    if (!input_signal || !filter_taps || !output_signal) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize input signal and filter taps with random values between 0.0 and 1.0
    for (int i = 0; i < signal_length; i++) {
        input_signal[i] = (double)mt_rand() / (double)UINT32_MAX;
    }

    for (int i = 0; i < num_filter_taps; i++) {
        filter_taps[i] = (double)mt_rand() / (double)UINT32_MAX;
    }
}

void run_computation() {
    // Apply the FIR filter (direct form convolution)
    for (int n = 0; n < signal_length; n++) {
        double sum = 0.0;
        for (int k = 0; k < num_filter_taps; k++) {
            // Convolve by multiplying filter tap with corresponding past input signal sample
            // This implicitly handles boundary conditions by only adding if the signal index is valid (>=0)
            if (n - k >= 0) {
                sum += filter_taps[k] * input_signal[n - k];
            }
        }
        output_signal[n] = sum;
    }

    // Accumulate the output to prevent dead code elimination
    accumulated_result = 0.0;
    for (int i = 0; i < signal_length; i++) {
        accumulated_result += output_signal[i];
    }
}

void cleanup() {
    free(input_signal);
    free(filter_taps);
    free(output_signal);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%f\n", accumulated_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}