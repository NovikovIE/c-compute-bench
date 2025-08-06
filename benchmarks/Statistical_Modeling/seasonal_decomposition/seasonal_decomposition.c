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
// --- END MERSENNE TWISTER ---

// Benchmark parameters
int num_time_steps;
int period;

// Data arrays
double *time_series;
double *trend;
double *seasonal;
double *residual;

// Temporary arrays for computation
double* seasonal_sum;
int* seasonal_count;
double* avg_seasonal;

// Final result accumulator
double final_result;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_time_steps> <period> <seed>\n", argv[0]);
        exit(1);
    }

    num_time_steps = atoi(argv[1]);
    period = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (num_time_steps <= 0 || period <= 1) {
        fprintf(stderr, "FATAL: num_time_steps and period must be positive (period > 1).\n");
        exit(1);
    }
    if (period >= num_time_steps / 2) {
        fprintf(stderr, "FATAL: period must be smaller than half of num_time_steps.\n");
        exit(1);
    }

    mt_seed(seed);

    time_series = (double *)malloc(num_time_steps * sizeof(double));
    trend = (double *)malloc(num_time_steps * sizeof(double));
    seasonal = (double *)malloc(num_time_steps * sizeof(double));
    residual = (double *)malloc(num_time_steps * sizeof(double));
    seasonal_sum = (double*)malloc(period * sizeof(double));
    seasonal_count = (int*)malloc(period * sizeof(int));
    avg_seasonal = (double*)malloc(period * sizeof(double));


    if (!time_series || !trend || !seasonal || !residual || !seasonal_sum || !seasonal_count || !avg_seasonal) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Generate synthetic time series data with trend, seasonality, and noise
    for (int i = 0; i < num_time_steps; i++) {
        double trend_component = 0.05 * i;
        double seasonal_component = 20.0 * sin(2.0 * M_PI * i / period);
        double noise_component = ((double)mt_rand() / (double)UINT32_MAX - 0.5) * 10.0;
        time_series[i] = trend_component + seasonal_component + noise_component;
    }
    
    final_result = 0.0;
}

void cleanup() {
    free(time_series);
    free(trend);
    free(seasonal);
    free(residual);
    free(seasonal_sum);
    free(seasonal_count);
    free(avg_seasonal);
}

void run_computation() {
    // ---- Step 1: Estimate Trend using a Centered Moving Average (naive implementation for computational intensity) ----
    int window_size = period;
    if (window_size % 2 == 0) {
        window_size++; // Ensure odd window size for simple centering
    }
    int half_window = window_size / 2;

    for (int i = half_window; i < num_time_steps - half_window; i++) {
        double sum = 0.0;
        for (int j = -half_window; j <= half_window; j++) {
            sum += time_series[i + j];
        }
        trend[i] = sum / window_size;
    }

    // Extrapolate trend to the edges by copying the first/last calculated values
    for (int i = 0; i < half_window; i++) {
        trend[i] = trend[half_window];
    }
    for (int i = num_time_steps - half_window; i < num_time_steps; i++) {
        trend[i] = trend[num_time_steps - half_window - 1];
    }
    
    // ---- Step 2: Estimate Seasonal Component ----
    for(int i = 0; i < period; ++i) {
        seasonal_sum[i] = 0.0;
        seasonal_count[i] = 0;
    }

    // Calculate detrended series and sum up for each seasonal period
    for (int i = 0; i < num_time_steps; i++) {
        double detrended_value = time_series[i] - trend[i];
        int p_idx = i % period;
        seasonal_sum[p_idx] += detrended_value;
        seasonal_count[p_idx]++;
    }

    // Calculate average seasonal effect
    double overall_seasonal_avg = 0.0;
    for (int i = 0; i < period; i++) {
        avg_seasonal[i] = (seasonal_count[i] > 0) ? seasonal_sum[i] / seasonal_count[i] : 0.0;
        overall_seasonal_avg += avg_seasonal[i];
    }
    overall_seasonal_avg /= period;

    // Normalize seasonal effects to sum to zero
    for (int i = 0; i < period; i++) {
        avg_seasonal[i] -= overall_seasonal_avg;
    }
    
    // ---- Step 3: Calculate final seasonal and residual series ----
    for (int i = 0; i < num_time_steps; i++) {
        seasonal[i] = avg_seasonal[i % period];
        residual[i] = time_series[i] - trend[i] - seasonal[i];
    }

    // ---- Step 4: Calculate final result (sum of absolute residuals) to prevent dead code elimination ----
    for (int i = 0; i < num_time_steps; i++) {
        final_result += fabs(residual[i]);
    }
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    // Print the final result to stdout to be checked
    printf("%f\n", final_result);

    // Calculate and print the execution time to stderr
    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
