#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <string.h>

// --- START: Mersenne Twister (MT19937) --- Do Not Modify ---
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
// --- END: Mersenne Twister (MT19937) ---

// Global data structure to hold all benchmark data
typedef struct {
    // Input parameter
    int num_tree_steps;

    // Financial model parameters
    double S0;   // Initial asset price
    double K;    // Strike price
    double T;    // Time to maturity (years)
    double r;    // Risk-free interest rate
    double sigma;// Volatility

    // Trinomial tree derived parameters
    double dt;   // Time step size
    double dx;   // Log-price step size
    double pu;   // Probability of up move
    double pd;   // Probability of down move
    double pm;   // Probability of middle move

    // Data arrays for computation
    double* V1; // Option values array 1
    double* V2; // Option values array 2

    // Final result
    double option_price;
} TrinomialBenchData;

TrinomialBenchData g_data;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_tree_steps> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_tree_steps = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);
    mt_seed(seed);
    
    if (g_data.num_tree_steps <= 0) {
        fprintf(stderr, "FATAL: num_tree_steps must be positive.\n");
        exit(1);
    }

    // Set fixed financial parameters for a European Call Option
    g_data.S0 = 100.0;    // Initial stock price
    g_data.K = 100.0;     // Strike price
    g_data.T = 1.0;       // Time to maturity in years
    g_data.r = 0.05;      // Risk-free rate
    g_data.sigma = 0.2;   // Volatility

    // Calculate derived trinomial tree parameters (Hull-White model)
    g_data.dt = g_data.T / g_data.num_tree_steps;
    g_data.dx = g_data.sigma * sqrt(3.0 * g_data.dt);
    double nu = g_data.r - 0.5 * g_data.sigma * g_data.sigma;

    double dx2 = g_data.dx * g_data.dx;
    g_data.pu = 0.5 * (((g_data.sigma * g_data.sigma * g_data.dt + nu * nu * g_data.dt * g_data.dt) / dx2) + (nu * g_data.dt / g_data.dx));
    g_data.pd = 0.5 * (((g_data.sigma * g_data.sigma * g_data.dt + nu * nu * g_data.dt * g_data.dt) / dx2) - (nu * g_data.dt / g_data.dx));
    g_data.pm = 1.0 - g_data.pu - g_data.pd;

    // Allocate memory for option value arrays
    // The number of nodes at the final step is 2*N+1
    int array_size = 2 * g_data.num_tree_steps + 1;
    g_data.V1 = (double*)malloc(array_size * sizeof(double));
    g_data.V2 = (double*)malloc(array_size * sizeof(double));

    if (!g_data.V1 || !g_data.V2) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }
}

void run_computation() {
    int N = g_data.num_tree_steps;
    double S0 = g_data.S0;
    double K = g_data.K;
    double dx = g_data.dx;
    double pu = g_data.pu;
    double pm = g_data.pm;
    double pd = g_data.pd;
    double discount = exp(-g_data.r * g_data.dt);

    // Use pointers for swapping arrays between steps efficiently
    double* P_next = g_data.V1;
    double* P_curr = g_data.V2;

    // 1. Calculate option values at maturity (time T, step N)
    // The states are k = -N, -N+1, ..., N. We map state k to array index j = k+N.
    for (int j = 0; j <= 2 * N; ++j) {
        int k = j - N; // State index
        double ST = S0 * exp(k * dx);
        P_next[j] = fmax(0.0, ST - K); // Call option payoff
    }

    // 2. Backward induction through the tree
    for (int i = N - 1; i >= 0; --i) {
        // At step i, states are k = -i, ..., i. Array index j = k+i.
        for (int j = 0; j <= 2 * i; ++j) {
            // The node (i, k) maps to P_curr[j]
            // It depends on nodes (i+1, k+1), (i+1, k), (i+1, k-1) from P_next
            // which map to array indices j+2, j+1, j
            P_curr[j] = discount * (pu * P_next[j + 2] + pm * P_next[j + 1] + pd * P_next[j]);
        }

        // Swap pointers for the next iteration
        double* temp = P_curr;
        P_curr = P_next;
        P_next = temp;
    }

    // The final result is at step 0, state 0 (array index 0)
    g_data.option_price = P_next[0];
}

void cleanup() {
    free(g_data.V1);
    free(g_data.V2);
    g_data.V1 = NULL;
    g_data.V2 = NULL;
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
    printf("%f\n", g_data.option_price);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
