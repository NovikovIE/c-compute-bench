#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- Mersenne Twister (MT19937) Generator (Verbatim) ---
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
// --- End of MT19937 ---

// Benchmark-specific global variables
int num_data_points;
double *data;
double *jackknife_estimates;
double final_result;

// Setup: parse arguments, allocate memory, and generate data
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_data_points> <seed>\n", argv[0]);
        exit(1);
    }

    num_data_points = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if (num_data_points <= 1) {
        fprintf(stderr, "Error: num_data_points must be greater than 1.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate memory
    data = (double *)malloc(num_data_points * sizeof(double));
    jackknife_estimates = (double *)malloc(num_data_points * sizeof(double));
    if (data == NULL || jackknife_estimates == NULL) {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        exit(1);
    }

    // Generate random data points
    for (int i = 0; i < num_data_points; ++i) {
        data[i] = (double)mt_rand() / (double)UINT32_MAX;
    }
}

// Core computation: Jackknife resampling to estimate the variance of the sample variance
void run_computation() {
    // The Jackknife procedure involves creating N 'leave-one-out' samples
    // and calculating a statistic for each. Here, the statistic is the sample variance.
    #pragma omp parallel for
    for (int i = 0; i < num_data_points; ++i) {
        // Calculate variance for the sample excluding point i
        double sum_x = 0.0;
        double sum_x_sq = 0.0;
        int n_minus_1 = num_data_points - 1;

        for (int j = 0; j < num_data_points; ++j) {
            if (i == j) continue;
            sum_x += data[j];
            sum_x_sq += data[j] * data[j];
        }

        double mean = sum_x / n_minus_1;
        // Sample variance denominator is (n-1) - 1, since the sample size is n-1
        double variance = (sum_x_sq / (n_minus_1 - 1.0)) - (mean * mean * n_minus_1) / (n_minus_1 - 1.0);
        jackknife_estimates[i] = variance;
    }

    // Calculate the mean of the N jackknife estimates
    double sum_estimates = 0.0;
    for (int i = 0; i < num_data_points; ++i) {
        sum_estimates += jackknife_estimates[i];
    }
    double mean_estimates = sum_estimates / num_data_points;

    // Calculate the sum of squared differences for the jackknife variance formula
    double sum_sq_diff_estimates = 0.0;
    for (int i = 0; i < num_data_points; ++i) {
        double diff = jackknife_estimates[i] - mean_estimates;
        sum_sq_diff_estimates += diff * diff;
    }

    // The final result is the Jackknife estimate of the variance of our statistic (sample variance)
    final_result = ((double)(num_data_points - 1) / num_data_points) * sum_sq_diff_estimates;
}

// Cleanup: free all allocated memory
void cleanup() {
    free(data);
    free(jackknife_estimates);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final result to stdout
    printf("%f\n", final_result);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
