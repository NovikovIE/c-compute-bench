#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- BEGIN: Mersenne Twister (DO NOT MODIFY) ---
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
// --- END: Mersenne Twister ---

// Global struct for benchmark data
static struct {
    int num_assets;
    int num_historical_points;
    uint32_t seed;

    double *historical_returns; // Flattened 2D array [asset][point]
    double *mean_returns;
    double *covariance_matrix; // Flattened 2D array [row][col]

    double final_result;
} g_data;

// Function to generate a random double between -0.05 and +0.05
double random_return() {
    return ((double)mt_rand() / (double)UINT32_MAX - 0.5) * 0.1;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_assets> <num_historical_points> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_assets = atoi(argv[1]);
    g_data.num_historical_points = atoi(argv[2]);
    g_data.seed = (uint32_t)atoi(argv[3]);

    mt_seed(g_data.seed);

    size_t returns_size = (size_t)g_data.num_assets * g_data.num_historical_points;
    g_data.historical_returns = (double*)malloc(returns_size * sizeof(double));
    g_data.mean_returns = (double*)malloc(g_data.num_assets * sizeof(double));
    g_data.covariance_matrix = (double*)malloc((size_t)g_data.num_assets * g_data.num_assets * sizeof(double));

    if (!g_data.historical_returns || !g_data.mean_returns || !g_data.covariance_matrix) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    for (size_t i = 0; i < returns_size; ++i) {
        g_data.historical_returns[i] = random_return();
    }
}

void run_computation() {
    int num_assets = g_data.num_assets;
    int num_points = g_data.num_historical_points;

    // 1. Calculate mean return for each asset
    for (int i = 0; i < num_assets; ++i) {
        double sum = 0.0;
        for (int p = 0; p < num_points; ++p) {
            sum += g_data.historical_returns[i * num_points + p];
        }
        g_data.mean_returns[i] = sum / num_points;
    }

    // 2. Calculate the covariance matrix
    // Cov(i, j) = E[(R_i - E[R_i]) * (R_j - E[R_j])]
    for (int i = 0; i < num_assets; ++i) {
        for (int j = i; j < num_assets; ++j) { // Matrix is symmetric, compute upper triangle
            double cov_sum = 0.0;
            double mean_i = g_data.mean_returns[i];
            double mean_j = g_data.mean_returns[j];
            for (int p = 0; p < num_points; ++p) {
                double dev_i = g_data.historical_returns[i * num_points + p] - mean_i;
                double dev_j = g_data.historical_returns[j * num_points + p] - mean_j;
                cov_sum += dev_i * dev_j;
            }
            // Use sample covariance (N-1 denominator)
            double covariance = cov_sum / (num_points > 1 ? (num_points - 1) : 1);
            g_data.covariance_matrix[i * num_assets + j] = covariance;
            if (i != j) {
                g_data.covariance_matrix[j * num_assets + i] = covariance; // And fill lower triangle
            }
        }
    }

    // 3. Calculate a final result to prevent dead code elimination
    double result_sum = 0.0;
    size_t cov_matrix_size = (size_t)num_assets * num_assets;
    for (size_t i = 0; i < cov_matrix_size; ++i) {
        result_sum += g_data.covariance_matrix[i];
    }
    g_data.final_result = result_sum;
}

void cleanup() {
    free(g_data.historical_returns);
    free(g_data.mean_returns);
    free(g_data.covariance_matrix);
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
    printf("%.12f\n", g_data.final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
