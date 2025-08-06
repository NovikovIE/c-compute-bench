#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- START MERSENNE TWISTER (DO NOT MODIFY) ---
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

// --- BENCHMARK DATA AND FUNCTIONS ---

typedef struct {
    // Parameters
    int duration_seconds;
    int sample_rate_hz;
    int filter_order;
    int num_modulation_steps;
    int num_harmonics;
    long num_samples;
    
    // Data arrays
    double* source_signal;
    double* filtered_signal;
    
    // Result
    double result_accumulator;
} BenchmarkData;

static BenchmarkData g_data;

// Setup: Parses arguments, allocates memory, and generates the initial rich audio signal.
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 7) {
        fprintf(stderr, "Usage: %s duration_seconds sample_rate_hz filter_order num_modulation_steps num_harmonics seed\n", argv[0]);
        exit(1);
    }
    
    g_data.duration_seconds = atoi(argv[1]);
    g_data.sample_rate_hz = atoi(argv[2]);
    g_data.filter_order = atoi(argv[3]);
    g_data.num_modulation_steps = atoi(argv[4]);
    g_data.num_harmonics = atoi(argv[5]);
    uint32_t seed = (uint32_t)atoi(argv[6]);
    
    mt_seed(seed);
    
    g_data.num_samples = (long)g_data.duration_seconds * g_data.sample_rate_hz;
    
    g_data.source_signal = (double*)malloc(g_data.num_samples * sizeof(double));
    g_data.filtered_signal = (double*)malloc(g_data.num_samples * sizeof(double));
    
    if (!g_data.source_signal || !g_data.filtered_signal) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Generate a harmonically rich source signal (sum of sine waves).
    // This simulates a classic subtractive synthesis source like a sawtooth wave.
    double base_freq = 110.0; // A2 note
    double* phases = (double*)malloc(g_data.num_harmonics * sizeof(double));
    for(int k=0; k < g_data.num_harmonics; ++k) {
        phases[k] = (mt_rand() / (double)UINT32_MAX) * 2.0 * M_PI; // Random initial phase
    }

    for (long i = 0; i < g_data.num_samples; ++i) {
        double sample = 0.0;
        double time = (double)i / g_data.sample_rate_hz;
        for (int k = 1; k <= g_data.num_harmonics; ++k) {
            sample += sin(2.0 * M_PI * base_freq * k * time + phases[k-1]);
        }
        g_data.source_signal[i] = sample / g_data.num_harmonics; // Normalize
        g_data.filtered_signal[i] = 0.0; // Initialize filtered signal
    }
    
    free(phases);
}

// Computation: Performs the core subtractive synthesis by applying a modulated filter.
void run_computation() {
    long samples_per_step = g_data.num_samples / g_data.num_modulation_steps;
    if (samples_per_step == 0) samples_per_step = 1;

    double* b = (double*)malloc((g_data.filter_order + 1) * sizeof(double)); // Feedforward coeffs
    double* a = (double*)malloc((g_data.filter_order + 1) * sizeof(double)); // Feedback coeffs
    if (!b || !a) {
        fprintf(stderr, "FATAL: Memory allocation failed in computation.\n");
        exit(1);
    }

    long current_sample_idx = 0;
    for (int step = 0; step < g_data.num_modulation_steps; ++step) {
        // "Recalculate" filter coefficients for this step. For this benchmark,
        // we generate random coefficients to simulate the workload of dynamically
        // adjusting a complex filter without the complexity of real filter design math.
        double a_sum = 0.0;
        for (int i = 0; i <= g_data.filter_order; ++i) {
            b[i] = (mt_rand() / (double)UINT32_MAX) * 2.0 - 1.0; // [-1.0, 1.0]
             // Keep 'a' coefficients small and positive to improve stability
            if (i > 0) {
                a[i] = (mt_rand() / (double)UINT32_MAX);
                a_sum += a[i];
            }
        }
        a[0] = 1.0; 

        // A crude attempt to prevent the filter from exploding by normalizing feedback terms.
        if (a_sum >= 1.0) {
            for (int i = 1; i <= g_data.filter_order; ++i) {
                a[i] /= (a_sum * 1.1);
            }
        }

        long end_sample = current_sample_idx + samples_per_step;
        if (step == g_data.num_modulation_steps - 1) {
            end_sample = g_data.num_samples; // Ensure all samples are processed
        }

        // Apply the IIR filter (Direct Form I) for this block of samples
        for (long i = current_sample_idx; i < end_sample; ++i) {
            double y_n = 0.0;
            // Feedforward part (from source signal)
            for (int j = 0; j <= g_data.filter_order; ++j) {
                if (i >= j) {
                    y_n += b[j] * g_data.source_signal[i - j];
                }
            }
            // Feedback part (from past output)
            for (int j = 1; j <= g_data.filter_order; ++j) {
                if (i >= j) {
                    y_n -= a[j] * g_data.filtered_signal[i - j];
                }
            }
            g_data.filtered_signal[i] = y_n; // No division as a[0] is 1.0
        }
        current_sample_idx = end_sample;
    }

    // Accumulate result to prevent dead code elimination.
    double sum = 0.0;
    for (long i = 0; i < g_data.num_samples; ++i) {
        sum += g_data.filtered_signal[i];
    }
    g_data.result_accumulator = sum;

    free(b);
    free(a);
}

// Cleanup: Frees memory allocated in setup.
void cleanup() {
    free(g_data.source_signal);
    free(g_data.filtered_signal);
}

// Main: Orchestrates the benchmark execution and timing.
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final accumulated result to stdout.
    printf("%.2f\n", g_data.result_accumulator);

    // Print the timing information to stderr.
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
