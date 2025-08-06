#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

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

// --- Benchmark Specific Code ---

typedef struct {
    int num_assets;
    int num_portfolios;
    double *asset_returns;
    double *covariance_matrix;
    double *portfolio_weights;
    double *portfolio_returns;
    double *portfolio_risks;
    double final_result;
} BenchmarkData;

static BenchmarkData g_data;

static double rand_double(double min, double max) {
    return min + ((double)mt_rand() / (double)UINT32_MAX) * (max - min);
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_assets> <num_portfolios> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_assets = atoi(argv[1]);
    g_data.num_portfolios = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);
    mt_seed(seed);

    int na = g_data.num_assets;
    int np = g_data.num_portfolios;

    g_data.asset_returns = (double *)malloc(na * sizeof(double));
    g_data.covariance_matrix = (double *)malloc(na * na * sizeof(double));
    g_data.portfolio_weights = (double *)malloc(np * na * sizeof(double));
    g_data.portfolio_returns = (double *)malloc(np * sizeof(double));
    g_data.portfolio_risks = (double *)malloc(np * sizeof(double));

    if (!g_data.asset_returns || !g_data.covariance_matrix || !g_data.portfolio_weights || !g_data.portfolio_returns || !g_data.portfolio_risks) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Generate asset characteristics
    double *asset_volatilities = (double *)malloc(na * sizeof(double));
    for (int i = 0; i < na; i++) {
        g_data.asset_returns[i] = rand_double(0.01, 0.15); // 1% to 15% annual return
        asset_volatilities[i] = rand_double(0.10, 0.40); // 10% to 40% annual volatility
    }

    // Generate a simplified covariance matrix
    for (int i = 0; i < na; i++) {
        for (int j = 0; j < na; j++) {
            if (i == j) {
                g_data.covariance_matrix[i * na + j] = asset_volatilities[i] * asset_volatilities[i];
            } else {
                double correlation = rand_double(-0.2, 0.7);
                g_data.covariance_matrix[i * na + j] = correlation * asset_volatilities[i] * asset_volatilities[j];
            }
        }
    }
    free(asset_volatilities);

    // Generate random portfolio weights
    for (int p = 0; p < np; p++) {
        double sum_of_weights = 0.0;
        for (int a = 0; a < na; a++) {
            double weight = rand_double(0.0, 1.0);
            g_data.portfolio_weights[p * na + a] = weight;
            sum_of_weights += weight;
        }
        // Normalize weights to sum to 1.0
        for (int a = 0; a < na; a++) {
            g_data.portfolio_weights[p * na + a] /= sum_of_weights;
        }
    }
}

void run_computation() {
    int na = g_data.num_assets;
    int np = g_data.num_portfolios;
    double final_checksum = 0.0;

    for (int p = 0; p < np; p++) {
        double portfolio_return = 0.0;
        double portfolio_variance = 0.0;

        // Calculate expected portfolio return: E(R_p) = W^T * E(R)
        for (int a = 0; a < na; a++) {
            portfolio_return += g_data.portfolio_weights[p * na + a] * g_data.asset_returns[a];
        }
        g_data.portfolio_returns[p] = portfolio_return;

        // Calculate portfolio variance: Var(R_p) = W^T * Cov * W
        for (int i = 0; i < na; i++) {
            for (int j = 0; j < na; j++) {
                portfolio_variance += g_data.portfolio_weights[p * na + i] * 
                                      g_data.portfolio_weights[p * na + j] * 
                                      g_data.covariance_matrix[i * na + j];
            }
        }
        g_data.portfolio_risks[p] = sqrt(portfolio_variance);

        final_checksum += portfolio_return - g_data.portfolio_risks[p];
    }

    g_data.final_result = final_checksum;
}

void cleanup() {
    free(g_data.asset_returns);
    free(g_data.covariance_matrix);
    free(g_data.portfolio_weights);
    free(g_data.portfolio_returns);
    free(g_data.portfolio_risks);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout to prevent dead code elimination
    printf("%f\n", g_data.final_result);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
