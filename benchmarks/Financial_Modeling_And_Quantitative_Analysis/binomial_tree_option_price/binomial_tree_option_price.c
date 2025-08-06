#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>
#include <string.h>

// --- Mersenne Twister (MT19937) Generator (Verbatim as Required) ---
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

// --- Benchmark Data and Globals ---
typedef struct {
    int num_tree_steps;
    double S0;    // Initial stock price
    double K;     // Strike price
    double T;     // Time to maturity (years)
    double r;     // Risk-free interest rate
    double sigma; // Volatility
    double* option_values; // Array for current time step's option values
    double result_price;   // Final computed option price
} BenchmarkData;

BenchmarkData g_data;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_tree_steps> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_tree_steps = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if (g_data.num_tree_steps <= 0) {
        fprintf(stderr, "ERROR: num_tree_steps must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    // Hardcoded financial parameters for consistency
    g_data.S0 = 100.0;    // Initial stock price
    g_data.K = 105.0;     // Strike price (out-of-the-money)
    g_data.T = 1.0;       // 1 year to maturity
    g_data.r = 0.05;      // 5% risk-free rate
    g_data.sigma = 0.2;   // 20% volatility

    // Allocate memory for the option values at each node of a time step
    g_data.option_values = (double*)malloc((g_data.num_tree_steps + 1) * sizeof(double));
    if (g_data.option_values == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    g_data.result_price = 0.0;
}

void run_computation() {
    // Binomial Tree parameters
    double dt = g_data.T / g_data.num_tree_steps;
    double u = exp(g_data.sigma * sqrt(dt));
    double d = 1.0 / u;
    double p = (exp(g_data.r * dt) - d) / (u - d);
    double discount_factor = exp(-g_data.r * dt);

    // 1. Calculate option values at the final nodes (at maturity)
    // This is for a European Call Option: max(S_T - K, 0)
    double current_stock_price = g_data.S0 * pow(d, g_data.num_tree_steps);
    double ud_ratio = u / d;
    for (int i = 0; i <= g_data.num_tree_steps; i++) {
        g_data.option_values[i] = fmax(0.0, current_stock_price - g_data.K);
        current_stock_price *= ud_ratio;
    }

    // 2. Work backwards through the tree to find the price at time 0
    for (int t = g_data.num_tree_steps - 1; t >= 0; t--) {
        for (int i = 0; i <= t; i++) {
            g_data.option_values[i] = discount_factor * (p * g_data.option_values[i + 1] + (1.0 - p) * g_data.option_values[i]);
        }
    }

    // 3. The final result is the value at the root of the tree
    g_data.result_price = g_data.option_values[0];
}

void cleanup() {
    if (g_data.option_values != NULL) {
        free(g_data.option_values);
        g_data.option_values = NULL;
    }
}

// --- Main Function ---
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%.10f\n", g_data.result_price);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
