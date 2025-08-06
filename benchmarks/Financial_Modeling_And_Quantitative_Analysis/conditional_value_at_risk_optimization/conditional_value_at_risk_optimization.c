/*
 * conditional_value_at_risk_optimization.c
 * 
 * A benchmark simulating a Conditional Value at Risk (CVaR) calculation for a financial portfolio.
 * CVaR (or Expected Shortfall) is a risk measure used to quantify the average loss
 * in the worst-case scenarios (e.g., the worst 5% of outcomes).
 *
 * The benchmark performs the following steps:
 * 1. Setup: Creates a large matrix of simulated asset returns for a given number
 *    of assets over a given number of scenarios. The returns are generated randomly.
 * 2. Computation: For a single, equally-weighted portfolio, it calculates the portfolio's
 *    loss for each scenario. It then sorts these losses to find the Value at Risk (VaR)
 *    threshold. Finally, it computes the CVaR by averaging all losses that exceed the VaR.
 * 3. Cleanup: Frees the memory used for the asset returns and portfolio losses.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- Mersenne Twister (MT19937) Generator - DO NOT MODIFY ---
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
// --- End of Mersenne Twister ---

// --- Benchmark Data and Globals ---

// Confidence level for VaR/CVaR calculation (e.g., 95%)
#define CONFIDENCE_LEVEL 0.95

typedef struct {
    size_t num_assets;
    size_t num_scenarios;
    double *returns; // Flattened 2D array [num_scenarios][num_assets]
    double *portfolio_losses;
    double final_cvar;
} BenchmarkData;

BenchmarkData *g_data = NULL;

// --- Utility Functions ---

// Comparison function for qsort
int compare_doubles(const void *a, const void *b) {
    double da = *(const double *)a;
    double db = *(const double *)b;
    if (da < db) return -1;
    if (da > db) return 1;
    return 0;
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_assets> <num_scenarios> <seed>\n", argv[0]);
        exit(1);
    }

    g_data = (BenchmarkData *)malloc(sizeof(BenchmarkData));
    if (!g_data) {
        perror("Failed to allocate memory for benchmark data");
        exit(1);
    }

    g_data->num_assets = strtol(argv[1], NULL, 10);
    g_data->num_scenarios = strtol(argv[2], NULL, 10);
    uint32_t seed = (uint32_t)strtol(argv[3], NULL, 10);

    mt_seed(seed);

    size_t total_returns = g_data->num_assets * g_data->num_scenarios;
    g_data->returns = (double *)malloc(total_returns * sizeof(double));
    g_data->portfolio_losses = (double *)malloc(g_data->num_scenarios * sizeof(double));

    if (!g_data->returns || !g_data->portfolio_losses) {
        perror("Failed to allocate memory for financial data");
        free(g_data->returns);
        free(g_data->portfolio_losses);
        free(g_data);
        exit(1);
    }

    // Generate plausible random asset returns (e.g., between -10% and +10%)
    for (size_t i = 0; i < total_returns; ++i) {
        double random_val = (double)mt_rand() / (double)UINT32_MAX;
        g_data->returns[i] = (random_val - 0.5) * 0.2;
    }
}

void run_computation() {
    const size_t num_assets = g_data->num_assets;
    const size_t num_scenarios = g_data->num_scenarios;
    const double equal_weight = 1.0 / (double)num_assets;

    // 1. Calculate the loss for each scenario for an equally-weighted portfolio
    for (size_t s = 0; s < num_scenarios; ++s) {
        double portfolio_return = 0.0;
        size_t scenario_offset = s * num_assets;
        for (size_t a = 0; a < num_assets; ++a) {
            portfolio_return += g_data->returns[scenario_offset + a] * equal_weight;
        }
        // Loss is the negative of the return
        g_data->portfolio_losses[s] = -portfolio_return;
    }

    // 2. Sort the losses to find the VaR threshold
    qsort(g_data->portfolio_losses, num_scenarios, sizeof(double), compare_doubles);

    // 3. Calculate CVaR: the average of the losses beyond the VaR threshold
    size_t var_index = (size_t)(num_scenarios * CONFIDENCE_LEVEL);
    if (var_index >= num_scenarios) {
        // Handle edge case where confidence level is too high or no scenarios exist.
        // The CVaR would be the single worst loss in this case.
        g_data->final_cvar = g_data->portfolio_losses[num_scenarios - 1];
        return;
    }

    double cvar_sum = 0.0;
    size_t tail_scenarios = num_scenarios - var_index;
    for (size_t i = var_index; i < num_scenarios; ++i) {
        cvar_sum += g_data->portfolio_losses[i];
    }

    g_data->final_cvar = cvar_sum / (double)tail_scenarios;
}

void cleanup() {
    if (g_data) {
        free(g_data->returns);
        free(g_data->portfolio_losses);
        free(g_data);
        g_data = NULL;
    }
}

// --- Main Execution --- 

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final result to stdout
    printf("%f\n", g_data->final_cvar);

    // Print the timing information to stderr
    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}