#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- MERSENNE TWISTER (DO NOT MODIFY) ---
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

// Global structure to hold benchmark data
typedef struct {
    long signal_length;
    int interpolation_factor;
    int filter_order;
    long upsampled_length;

    double* input_signal;
    double* upsampled_signal;
    double* output_signal;
    double* filter_coeffs;

    double result_checksum;
} BenchmarkData;

BenchmarkData g_data;

// Setup: parse arguments, allocate memory, and initialize data
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s signal_length interpolation_factor filter_order seed\n", argv[0]);
        exit(1);
    }

    g_data.signal_length = atol(argv[1]);
    g_data.interpolation_factor = atoi(argv[2]);
    g_data.filter_order = atoi(argv[3]);
    unsigned int seed = atoi(argv[4]);
    mt_seed(seed);

    g_data.upsampled_length = g_data.signal_length * g_data.interpolation_factor;

    // Allocate memory
    g_data.input_signal = (double*)malloc(g_data.signal_length * sizeof(double));
    g_data.upsampled_signal = (double*)malloc(g_data.upsampled_length * sizeof(double));
    g_data.output_signal = (double*)malloc(g_data.upsampled_length * sizeof(double));
    g_data.filter_coeffs = (double*)malloc(g_data.filter_order * sizeof(double));

    if (!g_data.input_signal || !g_data.upsampled_signal || !g_data.output_signal || !g_data.filter_coeffs) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize input signal with random values between -1.0 and 1.0
    for (long i = 0; i < g_data.signal_length; i++) {
        g_data.input_signal[i] = ((double)mt_rand() / (double)UINT32_MAX) * 2.0 - 1.0;
    }

    // Initialize a simple low-pass filter (moving average)
    double coeff_val = 1.0 / g_data.filter_order;
    for (int i = 0; i < g_data.filter_order; i++) {
        g_data.filter_coeffs[i] = coeff_val;
    }

    g_data.result_checksum = 0.0;
}

// Computation: perform upsampling and filtering
void run_computation() {
    // Step 1: Upsample by inserting zeros (zero-stuffing)
    for (long i = 0; i < g_data.signal_length; ++i) {
        g_data.upsampled_signal[i * g_data.interpolation_factor] = g_data.input_signal[i];
        for (int j = 1; j < g_data.interpolation_factor; ++j) {
            g_data.upsampled_signal[i * g_data.interpolation_factor + j] = 0.0;
        }
    }

    // Step 2: Apply the FIR filter (convolution) for interpolation
    for (long i = 0; i < g_data.upsampled_length; ++i) {
        double sum = 0.0;
        for (int j = 0; j < g_data.filter_order; ++j) {
            long k = i - j;
            if (k >= 0) {
                sum += g_data.upsampled_signal[k] * g_data.filter_coeffs[j];
            }
        }
        g_data.output_signal[i] = sum;
    }

    // Step 3: Compute a checksum to prevent dead-code elimination
    double checksum = 0.0;
    for (long i = 0; i < g_data.upsampled_length; ++i) {
        checksum += g_data.output_signal[i];
    }
    g_data.result_checksum = checksum;
}

// Cleanup: free all allocated memory
void cleanup() {
    free(g_data.input_signal);
    free(g_data.upsampled_signal);
    free(g_data.output_signal);
    free(g_data.filter_coeffs);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final checksum to stdout
    printf("%f\n", g_data.result_checksum);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
