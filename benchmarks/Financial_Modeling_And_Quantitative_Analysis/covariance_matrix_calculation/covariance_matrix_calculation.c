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

// --- Benchmark Globals ---
static int num_assets;
static int num_historical_points;
static double* asset_returns;       // 2D Array: num_assets x num_historical_points
static double* mean_returns;        // 1D Array: num_assets
static double* covariance_matrix;   // 2D Array: num_assets x num_assets
static double result_accumulator = 0.0;

// Macros for 2D array access on 1D heap allocation
#define RETURNS(asset_idx, point_idx) asset_returns[(asset_idx) * num_historical_points + (point_idx)]
#define COV_MATRIX(row, col) covariance_matrix[(row) * num_assets + (col)]

// Function to generate a random double between -0.1 and 0.1
double random_return() {
    // Generates a value in [0, 1] then scales and shifts it to [-0.1, 0.1]
    return ((double)mt_rand() / (double)UINT32_MAX) * 0.2 - 0.1;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_assets> <num_historical_points> <seed>\n", argv[0]);
        exit(1);
    }

    num_assets = atoi(argv[1]);
    num_historical_points = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    mt_seed(seed);

    // Allocate memory on the heap
    asset_returns = (double*)malloc((size_t)num_assets * num_historical_points * sizeof(double));
    mean_returns = (double*)malloc((size_t)num_assets * sizeof(double));
    covariance_matrix = (double*)malloc((size_t)num_assets * num_assets * sizeof(double));

    if (!asset_returns || !mean_returns || !covariance_matrix) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Initialize asset returns with random data
    for (int i = 0; i < num_assets; ++i) {
        for (int j = 0; j < num_historical_points; ++j) {
            RETURNS(i, j) = random_return();
        }
    }
}

void run_computation() {
    // Step 1: Calculate the mean return for each asset
    for (int i = 0; i < num_assets; ++i) {
        double sum = 0.0;
        for (int j = 0; j < num_historical_points; ++j) {
            sum += RETURNS(i, j);
        }
        mean_returns[i] = sum / num_historical_points;
    }

    // Step 2: Calculate the covariance matrix
    // The covariance of asset i and asset j is:
    // Cov(i, j) = sum[(R_i,t - E[R_i]) * (R_j,t - E[R_j])] / (N - 1)
    // where t is a historical point and N is num_historical_points.
    for (int i = 0; i < num_assets; ++i) {
        // The matrix is symmetric, so we only need to compute the upper triangle
        for (int j = i; j < num_assets; ++j) {
            double cov_sum = 0.0;
            for (int k = 0; k < num_historical_points; ++k) {
                double deviation_i = RETURNS(i, k) - mean_returns[i];
                double deviation_j = RETURNS(j, k) - mean_returns[j];
                cov_sum += deviation_i * deviation_j;
            }
            double covariance = cov_sum / (num_historical_points - 1);
            COV_MATRIX(i, j) = covariance;
            COV_MATRIX(j, i) = covariance; // Assign to the lower triangle as well
        }
    }

    // Step 3: Accumulate a result to prevent dead code elimination
    for (int i = 0; i < num_assets; ++i) {
        for (int j = 0; j < num_assets; ++j) {
            result_accumulator += COV_MATRIX(i, j);
        }
    }
}

void cleanup() {
    free(asset_returns);
    free(mean_returns);
    free(covariance_matrix);
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
    printf("%f\n", result_accumulator);

    // Print timing info to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
