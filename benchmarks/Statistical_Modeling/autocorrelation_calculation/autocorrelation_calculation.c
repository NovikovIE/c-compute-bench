/*
 * BENCHMARK: Statistical Modeling - Autocorrelation Calculation
 * This benchmark computes the autocorrelation function for a randomly generated
 * time series. Autocorrelation is the correlation of a signal with a delayed
 * copy of itself, used frequently in time series analysis and signal processing.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// MERSENNE TWISTER (MT19937) GENERATOR - DO NOT MODIFY
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
// END OF MT19937

// --- BENCHMARK DATA AND PARAMETERS ---
static int num_time_steps;
static int max_lag;

static double *time_series;
static double *autocorrelation_results;

// Final result accumulator to prevent dead code elimination
static double final_result_sum;

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_time_steps> <max_lag> <seed>\n", argv[0]);
        exit(1);
    }

    num_time_steps = atoi(argv[1]);
    max_lag = atoi(argv[2]);
    uint32_t seed = (uint32_t)strtoul(argv[3], NULL, 10);

    if (num_time_steps <= 0 || max_lag < 0 || max_lag >= num_time_steps) {
        fprintf(stderr, "Invalid arguments: num_time_steps must be > 0 and max_lag must be >= 0 and < num_time_steps.\n");
        exit(1);
    }

    mt_seed(seed);

    time_series = (double *)malloc(num_time_steps * sizeof(double));
    if (time_series == NULL) {
        fprintf(stderr, "Failed to allocate memory for time_series.\n");
        exit(1);
    }

    // max_lag+1 to include lag 0
    autocorrelation_results = (double *)malloc((max_lag + 1) * sizeof(double));
    if (autocorrelation_results == NULL) {
        fprintf(stderr, "Failed to allocate memory for autocorrelation_results.\n");
        free(time_series);
        exit(1);
    }

    // Generate a random time series with values between 0.0 and 1.0
    for (int i = 0; i < num_time_steps; ++i) {
        time_series[i] = (double)mt_rand() / (double)UINT32_MAX;
    }
}

void run_computation() {
    // 1. Calculate the mean of the time series
    double mean = 0.0;
    for (int i = 0; i < num_time_steps; ++i) {
        mean += time_series[i];
    }
    mean /= (double)num_time_steps;

    // 2. Calculate the variance of the series (denominator of the autocorrelation formula)
    double denominator = 0.0;
    for (int i = 0; i < num_time_steps; ++i) {
        denominator += (time_series[i] - mean) * (time_series[i] - mean);
    }

    // Avoid division by zero if variance is zero (e.g., constant series)
    if (denominator == 0.0) {
        for (int k = 0; k <= max_lag; ++k) {
            autocorrelation_results[k] = 1.0;
        }
    } else {
        // 3. Calculate autocorrelation for each lag from 0 to max_lag
        for (int k = 0; k <= max_lag; ++k) {
            double numerator = 0.0;
            for (int t = 0; t < num_time_steps - k; ++t) {
                numerator += (time_series[t] - mean) * (time_series[t + k] - mean);
            }
            autocorrelation_results[k] = numerator / denominator;
        }
    }

    // 4. Sum up the results to produce a single final value
    final_result_sum = 0.0;
    for (int k = 0; k <= max_lag; ++k) {
        final_result_sum += autocorrelation_results[k];
    }
}

void cleanup() {
    free(time_series);
    free(autocorrelation_results);
}

// --- MAIN FUNCTION ---

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout
    printf("%f\n", final_result_sum);

    // Print timing information to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
