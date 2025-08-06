#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

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


// --- Benchmark Globals ---
typedef struct {
    int num_assets;
    int num_simulations;
    int time_horizon_days;

    double* initial_prices;      // Initial price for each asset
    double* drifts;              // Annualized drift (expected return)
    double* volatilities;        // Annualized volatility
    double* weights;             // Number of shares for each asset

    double* simulation_pnl;      // Profit/Loss for each simulation path
    double final_var_result;     // Final calculated 99% Value at Risk
} BenchmarkData;

BenchmarkData g_data;

// --- Utility Functions ---

// Generates a random double in [0, 1) from the MT19937 generator
double rand_double() {
    return (double)mt_rand() / 4294967296.0;
}

// Generates a standard normal random variable using the Box-Muller transform
double generate_normal() {
    static int has_spare = 0;
    static double spare;

    if (has_spare) {
        has_spare = 0;
        return spare;
    }

    has_spare = 1;
    double u1, u2;
    do {
        u1 = rand_double();
        u2 = rand_double();
    } while (u1 <= 1e-9); // Avoid log(0)

    double r = sqrt(-2.0 * log(u1));
    double theta = 2.0 * M_PI * u2;

    spare = r * sin(theta);
    return r * cos(theta);
}

// Comparison function for qsort
int compare_doubles(const void* a, const void* b) {
    double da = *(const double*)a;
    double db = *(const double*)b;
    if (da < db) return -1;
    if (da > db) return 1;
    return 0;
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_assets> <num_simulations> <time_horizon_days> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_assets = atoi(argv[1]);
    g_data.num_simulations = atoi(argv[2]);
    g_data.time_horizon_days = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    mt_seed(seed);

    // Allocate memory
    g_data.initial_prices = (double*)malloc(g_data.num_assets * sizeof(double));
    g_data.drifts = (double*)malloc(g_data.num_assets * sizeof(double));
    g_data.volatilities = (double*)malloc(g_data.num_assets * sizeof(double));
    g_data.weights = (double*)malloc(g_data.num_assets * sizeof(double));
    g_data.simulation_pnl = (double*)malloc(g_data.num_simulations * sizeof(double));

    if (!g_data.initial_prices || !g_data.drifts || !g_data.volatilities || !g_data.weights || !g_data.simulation_pnl) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        exit(1);
    }

    // Initialize data with random values
    for (int i = 0; i < g_data.num_assets; ++i) {
        g_data.initial_prices[i] = 50.0 + 100.0 * rand_double(); // Price between 50 and 150
        g_data.drifts[i] = 0.01 + 0.14 * rand_double();          // Annual drift between 1% and 15%
        g_data.volatilities[i] = 0.10 + 0.30 * rand_double();     // Annual volatility between 10% and 40%
        g_data.weights[i] = 1.0;                                  // Assume portfolio holds 1 share of each asset
    }
}

void run_computation() {
    const double dt = 1.0 / 252.0; // Time step for one day (assuming 252 trading days/year)
    const double sqrt_dt = sqrt(dt);
    double* current_prices = (double*)malloc(g_data.num_assets * sizeof(double)); 
    if (!current_prices) {
        fprintf(stderr, "Error: Could not allocate memory for current_prices\n");
        exit(1);
    }

    // Calculate initial portfolio value
    double initial_portfolio_value = 0.0;
    for (int i = 0; i < g_data.num_assets; ++i) {
        initial_portfolio_value += g_data.initial_prices[i] * g_data.weights[i];
    }

    // Main Monte Carlo simulation loop
    for (int s = 0; s < g_data.num_simulations; ++s) {
        // Reset current prices to initial prices for this simulation run
        for (int i = 0; i < g_data.num_assets; ++i) {
            current_prices[i] = g_data.initial_prices[i];
        }

        // Simulate daily price changes for the time horizon
        for (int d = 0; d < g_data.time_horizon_days; ++d) {
            for (int a = 0; a < g_data.num_assets; ++a) {
                double normal_random = generate_normal();
                double drift_term = (g_data.drifts[a] - 0.5 * g_data.volatilities[a] * g_data.volatilities[a]) * dt;
                double diffusion_term = g_data.volatilities[a] * sqrt_dt * normal_random;
                current_prices[a] *= exp(drift_term + diffusion_term);
            }
        }

        // Calculate final portfolio value for this simulation
        double final_portfolio_value = 0.0;
        for (int i = 0; i < g_data.num_assets; ++i) {
            final_portfolio_value += current_prices[i] * g_data.weights[i];
        }

        // Store the profit or loss (P/L)
        g_data.simulation_pnl[s] = final_portfolio_value - initial_portfolio_value;
    }

    // Sort the P/L results to find the percentile
    qsort(g_data.simulation_pnl, g_data.num_simulations, sizeof(double), compare_doubles);

    // Calculate Value at Risk (VaR) at 99% confidence level
    // This is the P/L at the 1st percentile of the loss distribution
    int var_index = (int)(g_data.num_simulations * 0.01);
    g_data.final_var_result = g_data.simulation_pnl[var_index];

    free(current_prices);
}

void cleanup() {
    free(g_data.initial_prices);
    free(g_data.drifts);
    free(g_data.volatilities);
    free(g_data.weights);
    free(g_data.simulation_pnl);
}

int main(int argc, char* argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final result to stdout
    printf("%.4f\n", g_data.final_var_result); 

    // Print the execution time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
