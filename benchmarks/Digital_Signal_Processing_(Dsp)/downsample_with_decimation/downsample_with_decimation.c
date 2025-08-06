#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- Mersenne Twister (Do Not Modify) ---
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
// --- End Mersenne Twister ---

// Benchmark parameters and data
size_t g_signal_length;
int g_decimation_factor;
int g_filter_order;
double* g_input_signal;
double* g_filter_coeffs;
double* g_output_signal;
double g_final_result = 0.0;

// Function to generate a random double between 0.0 and 1.0
double rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s signal_length decimation_factor filter_order seed\n", argv[0]);
        exit(1);
    }

    g_signal_length = (size_t)atoll(argv[1]);
    g_decimation_factor = atoi(argv[2]);
    g_filter_order = atoi(argv[3]);
    uint32_t seed = (uint32_t)strtoul(argv[4], NULL, 10);

    if (g_decimation_factor <= 0 || g_filter_order <= 0 || g_signal_length <= 0) {
        fprintf(stderr, "ERROR: Invalid parameters.\n");
        exit(1);
    }
    if (g_filter_order > g_signal_length) {
        fprintf(stderr, "ERROR: Filter order cannot be greater than signal length.\n");
        exit(1);
    }

    mt_seed(seed);

    g_input_signal = (double*)malloc(g_signal_length * sizeof(double));
    g_filter_coeffs = (double*)malloc(g_filter_order * sizeof(double));
    size_t output_length = g_signal_length / g_decimation_factor;
    g_output_signal = (double*)malloc(output_length * sizeof(double));

    if (!g_input_signal || !g_filter_coeffs || !g_output_signal) {
        fprintf(stderr, "ERROR: Memory allocation failed.\n");
        exit(1);
    }

    // Generate a random input signal from -1.0 to 1.0
    for (size_t i = 0; i < g_signal_length; ++i) {
        g_input_signal[i] = rand_double() * 2.0 - 1.0; 
    }

    // Generate random filter coefficients for a basic low-pass filter and normalize them
    double coeff_sum = 0.0;
    for (int i = 0; i < g_filter_order; ++i) {
        g_filter_coeffs[i] = rand_double();
        coeff_sum += g_filter_coeffs[i];
    }
    // Normalize to make it a weighted average
    for (int i = 0; i < g_filter_order; ++i) {
        g_filter_coeffs[i] /= coeff_sum;
    }
}

void run_computation() {
    size_t output_length = g_signal_length / g_decimation_factor;
    double accumulator = 0.0;
    
    // Apply a causal FIR filter (convolution) and decimate in one pass
    for (size_t i = 0; i < output_length; ++i) {
        size_t input_index = i * g_decimation_factor;
        double sample_val = 0.0;
        
        // This inner loop performs the convolution for a single output sample.
        for (int j = 0; j < g_filter_order; ++j) {
            // Check boundary to avoid reading before the start of the signal array.
            if (input_index >= j) {
                sample_val += g_input_signal[input_index - j] * g_filter_coeffs[j];
            }
        }
        g_output_signal[i] = sample_val;
    }

    // Accumulate the result to prevent dead code elimination by the compiler
    for (size_t i = 0; i < output_length; i++) {
        accumulator += g_output_signal[i];
    }
    g_final_result = accumulator;
}

void cleanup() {
    free(g_input_signal);
    free(g_filter_coeffs);
    free(g_output_signal);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();
    
    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout
    printf("%f\n", g_final_result);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
