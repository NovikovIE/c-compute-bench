#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- Mersenne Twister Generator (Included Verbatim) ---
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

// --- Benchmark Globals ---

// Parameters
long g_signal_length;
int g_num_feedback_taps; // N from the difference equation
int g_num_feedforward_taps; // M from the difference equation

// Data arrays
double* g_input_signal = NULL;
double* g_output_signal = NULL;
double* g_feedback_coeffs = NULL; // 'a' coefficients
double* g_feedforward_coeffs = NULL; // 'b' coefficients

// Final result accumulator
double g_final_sum = 0.0;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <signal_length> <num_feedback_taps> <num_feedforward_taps> <seed>\n", argv[0]);
        exit(1);
    }

    g_signal_length = atol(argv[1]);
    g_num_feedback_taps = atoi(argv[2]);
    g_num_feedforward_taps = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    mt_seed(seed);

    g_input_signal = (double*)malloc(g_signal_length * sizeof(double));
    g_output_signal = (double*)calloc(g_signal_length, sizeof(double)); // calloc initializes to zero
    g_feedback_coeffs = (double*)malloc(g_num_feedback_taps * sizeof(double));
    g_feedforward_coeffs = (double*)malloc(g_num_feedforward_taps * sizeof(double));

    if (!g_input_signal || !g_output_signal || !g_feedback_coeffs || !g_feedforward_coeffs) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Generate random input signal
    for (long i = 0; i < g_signal_length; ++i) {
        g_input_signal[i] = (mt_rand() / (double)UINT32_MAX) * 2.0 - 1.0;
    }

    // Generate random feedforward coefficients
    for (int i = 0; i < g_num_feedforward_taps; ++i) {
        g_feedforward_coeffs[i] = (mt_rand() / (double)UINT32_MAX) * 2.0 - 1.0;
    }

    // Generate feedback coefficients. Must have a_0 = 1.
    // Keep other a_k small to increase chances of a stable filter.
    g_feedback_coeffs[0] = 1.0;
    for (int i = 1; i < g_num_feedback_taps; ++i) {
        g_feedback_coeffs[i] = (mt_rand() / (double)UINT32_MAX) * 0.1 - 0.05;
    }
}

void run_computation() {
    // IIR filter - Direct Form I implementation
    // y[n] = ( sum_{k=0 to M} b_k*x[n-k] - sum_{k=1 to N} a_k*y[n-k] ) / a_0
    // We assume a_0 = 1.

    for (long n = 0; n < g_signal_length; ++n) {
        double y_n = 0.0;

        // Feedforward part: sum_{k=0 to M} b_k*x[n-k]
        for (int k = 0; k < g_num_feedforward_taps; ++k) {
            if (n >= k) {
                y_n += g_feedforward_coeffs[k] * g_input_signal[n - k];
            }
        }

        // Feedback part: - sum_{k=1 to N} a_k*y[n-k]
        for (int k = 1; k < g_num_feedback_taps; ++k) {
            if (n >= k) {
                y_n -= g_feedback_coeffs[k] * g_output_signal[n - k];
            }
        }

        g_output_signal[n] = y_n;
    }

    // Accumulate result to prevent dead code elimination
    double sum = 0.0;
    for (long i = 0; i < g_signal_length; ++i) {
        sum += g_output_signal[i];
    }
    g_final_sum = sum;
}

void cleanup() {
    free(g_input_signal);
    free(g_output_signal);
    free(g_feedback_coeffs);
    free(g_feedforward_coeffs);
    g_input_signal = NULL;
    g_output_signal = NULL;
    g_feedback_coeffs = NULL;
    g_feedforward_coeffs = NULL;
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    // The result is printed to stdout
    printf("%f\n", g_final_sum);

    cleanup();

    // The time is printed to stderr
    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
