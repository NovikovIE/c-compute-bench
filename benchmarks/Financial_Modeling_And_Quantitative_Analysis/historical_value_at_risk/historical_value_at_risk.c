#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- Mersenne Twister (MT19937) Generator ---
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
// --- End of Mersenne Twister ---


// --- Benchmark Data Structures ---
#define CONFIDENCE_LEVEL 0.99

// A struct to hold all benchmark data.
typedef struct {
    int num_assets;
    int num_historical_days;

    // Input data
    double** historical_returns; // Matrix of [day][asset] returns
    double* portfolio_weights;   // Weight of each asset in the portfolio

    // Workspace and result
    double* portfolio_daily_returns; // Computed daily return of the whole portfolio
    double final_VaR;                // The final calculated Value-at-Risk
} benchmark_data_t;

// Global pointer to the benchmark data
benchmark_data_t* g_data = NULL;


// --- Benchmark Functions ---

// Comparison function for qsort
int compare_doubles(const void* a, const void* b) {
    double da = *(const double*)a;
    double db = *(const double*)b;
    if (da < db) return -1;
    if (da > db) return 1;
    return 0;
}

// Generate a random daily return, e.g., between -5% and +5%
double generate_random_return() {
    return ((double)mt_rand() / (double)UINT32_MAX - 0.5) * 0.1;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_assets> <num_historical_days> <seed>\n", argv[0]);
        exit(1);
    }

    g_data = (benchmark_data_t*)malloc(sizeof(benchmark_data_t));
    if (!g_data) {
        perror("Failed to allocate memory for benchmark_data_t");
        exit(1);
    }

    g_data->num_assets = atoi(argv[1]);
    g_data->num_historical_days = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);
    mt_seed(seed);

    // Allocate memory for historical returns matrix
    g_data->historical_returns = (double**)malloc(g_data->num_historical_days * sizeof(double*));
    if (!g_data->historical_returns) {
         perror("Failed to allocate memory for historical_returns rows");
         exit(1);
    }
    for (int i = 0; i < g_data->num_historical_days; ++i) {
        g_data->historical_returns[i] = (double*)malloc(g_data->num_assets * sizeof(double));
        if (!g_data->historical_returns[i]) {
            perror("Failed to allocate memory for historical_returns columns");
            exit(1);
        }
    }

    // Allocate other arrays
    g_data->portfolio_weights = (double*)malloc(g_data->num_assets * sizeof(double));
    g_data->portfolio_daily_returns = (double*)malloc(g_data->num_historical_days * sizeof(double));

    if (!g_data->portfolio_weights || !g_data->portfolio_daily_returns) {
        perror("Failed to allocate memory for portfolio arrays");
        exit(1);
    }

    // Populate the historical returns with random data
    for (int d = 0; d < g_data->num_historical_days; ++d) {
        for (int a = 0; a < g_data->num_assets; ++a) {
            g_data->historical_returns[d][a] = generate_random_return();
        }
    }

    // Populate portfolio weights (using an equally weighted portfolio for simplicity)
    double equal_weight = 1.0 / g_data->num_assets;
    for (int a = 0; a < g_data->num_assets; ++a) {
        g_data->portfolio_weights[a] = equal_weight;
    }
    
    g_data->final_VaR = 0.0;
}

void run_computation() {
    // Step 1: Calculate the historical daily returns of the entire portfolio.
    for (int d = 0; d < g_data->num_historical_days; ++d) {
        double daily_return = 0.0;
        for (int a = 0; a < g_data->num_assets; ++a) {
            daily_return += g_data->historical_returns[d][a] * g_data->portfolio_weights[a];
        }
        g_data->portfolio_daily_returns[d] = daily_return;
    }

    // Step 2: Sort the portfolio's daily returns in ascending order.
    // The worst losses will be at the beginning of the array.
    qsort(g_data->portfolio_daily_returns, g_data->num_historical_days, sizeof(double), compare_doubles);

    // Step 3: Determine the Value-at-Risk (VaR) at the specified confidence level.
    // For a 99% confidence level, we find the return at the 1st percentile.
    int var_index = (int)((1.0 - CONFIDENCE_LEVEL) * g_data->num_historical_days);

    // A check to ensure index is within bounds
    if (var_index < 0) var_index = 0;
    if (var_index >= g_data->num_historical_days) var_index = g_data->num_historical_days - 1;

    // VaR is conventionally reported as a positive number representing the potential loss,
    // so we negate the (typically negative) return value from the sorted list.
    g_data->final_VaR = -g_data->portfolio_daily_returns[var_index];
}

void cleanup() {
    if (!g_data) return;

    for (int i = 0; i < g_data->num_historical_days; ++i) {
        free(g_data->historical_returns[i]);
    }
    free(g_data->historical_returns);
    free(g_data->portfolio_weights);
    free(g_data->portfolio_daily_returns);
    free(g_data);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
    // Print the final result to stdout to be checked for correctness.
    printf("%f\n", g_data->final_VaR);

    cleanup();

    // Print the timing information to stderr for performance measurement.
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
