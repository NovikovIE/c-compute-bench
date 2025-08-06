#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <string.h>

// Define M_PI if not provided by math.h
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- Mersenne Twister (MT19937) Generator - Do Not Modify ---
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

// --- Global Data Structure ---
typedef struct {
    // Parameters from argv
    int num_simulations;
    int num_time_steps;
    int polynomial_degree;
    uint32_t seed;

    // Fixed financial model parameters
    double S0;   // Initial asset price
    double K;    // Strike price
    double r;    // Risk-free interest rate
    double sigma;// Volatility
    double T;    // Time to maturity (in years)

    // Data allocated on the heap
    double **asset_prices; // Matrix of size [num_simulations][num_time_steps + 1]
    double *cash_flows;    // Vector of size [num_simulations]

    // Final computed result
    double option_price;
} BenchmarkData;

BenchmarkData g_data;

// --- Helper Functions ---
double generate_normal_dist() {
    // Box-Muller transform to generate standard normal random variables.
    // Uses two uniform random numbers to generate two normal numbers.
    // Caches the second number for the next call.
    static int has_cached = 0;
    static double cached_val;
    if (has_cached) {
        has_cached = 0;
        return cached_val;
    }

    // Generate U1, U2 from (0, 1], not including 0 to avoid log(0)
    double u1 = (mt_rand() + 1.0) / (UINT32_MAX + 2.0);
    double u2 = (mt_rand() + 1.0) / (UINT32_MAX + 2.0);

    double log_u1 = sqrt(-2.0 * log(u1));
    double two_pi_u2 = 2.0 * M_PI * u2;

    cached_val = log_u1 * sin(two_pi_u2);
    has_cached = 1;
    
    return log_u1 * cos(two_pi_u2);
}

// Simple Gaussian elimination without pivoting to solve Ax = b
// (A is n x n, b is n x 1, x is n x 1)
int solve_linear_system(double **A, double *b, double *x, int n) {
    double **aug = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) {
        aug[i] = (double *)malloc((n + 1) * sizeof(double));
        for (int j = 0; j < n; j++) {
            aug[i][j] = A[i][j];
        }
        aug[i][n] = b[i];
    }

    // Forward elimination
    for (int i = 0; i < n; i++) {
        int max_row = i;
        for (int k = i + 1; k < n; k++) {
             if (fabs(aug[k][i]) > fabs(aug[max_row][i])) {
                max_row = k;
            }
        }
        double* temp = aug[i];
        aug[i] = aug[max_row];
        aug[max_row] = temp;

        if (fabs(aug[i][i]) < 1e-12) { // Singular or nearly singular
             for(int k=0; k<n; k++) free(aug[k]);
             free(aug);
             return -1; // Failure
        }
        for (int k = i + 1; k < n; k++) {
            double factor = aug[k][i] / aug[i][i];
            for (int j = i; j < n + 1; j++) {
                aug[k][j] -= factor * aug[i][j];
            }
        }
    }

    // Back substitution
    for (int i = n - 1; i >= 0; i--) {
        x[i] = aug[i][n];
        for (int j = i + 1; j < n; j++) {
            x[i] -= aug[i][j] * x[j];
        }
        x[i] /= aug[i][i];
    }
    
    for(int i=0; i<n; i++) free(aug[i]);
    free(aug);
    return 0; // Success
}

// --- Benchmark Functions ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_simulations> <num_time_steps> <polynomial_degree> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_simulations = atoi(argv[1]);
    g_data.num_time_steps = atoi(argv[2]);
    g_data.polynomial_degree = atoi(argv[3]);
    g_data.seed = (uint32_t)strtoul(argv[4], NULL, 10);

    mt_seed(g_data.seed);

    // Set fixed financial parameters
    g_data.S0 = 100.0;
    g_data.K = 100.0;
    g_data.r = 0.05;
    g_data.sigma = 0.2;
    g_data.T = 1.0;

    // Allocate memory
    g_data.asset_prices = (double **)malloc(g_data.num_simulations * sizeof(double *));
    for (int i = 0; i < g_data.num_simulations; i++) {
        g_data.asset_prices[i] = (double *)malloc((g_data.num_time_steps + 1) * sizeof(double));
    }
    g_data.cash_flows = (double *)malloc(g_data.num_simulations * sizeof(double));

    // Generate asset price paths using Geometric Brownian Motion
    double dt = g_data.T / g_data.num_time_steps;
    double drift = (g_data.r - 0.5 * g_data.sigma * g_data.sigma) * dt;
    double diffusion = g_data.sigma * sqrt(dt);

    for (int i = 0; i < g_data.num_simulations; i++) {
        g_data.asset_prices[i][0] = g_data.S0;
        for (int j = 0; j < g_data.num_time_steps; j++) {
            double Z = generate_normal_dist();
            g_data.asset_prices[i][j+1] = g_data.asset_prices[i][j] * exp(drift + diffusion * Z);
        }
    }
}

void run_computation() {
    double dt = g_data.T / g_data.num_time_steps;
    double discount_factor = exp(-g_data.r * dt);
    int d = g_data.polynomial_degree;

    // Initialize cash flows at maturity (T) for an American Put option
    for (int i = 0; i < g_data.num_simulations; i++) {
        double stock_price = g_data.asset_prices[i][g_data.num_time_steps];
        double payoff = (g_data.K - stock_price > 0) ? (g_data.K - stock_price) : 0.0;
        g_data.cash_flows[i] = payoff;
    }

    // Backward induction using Longstaff-Schwartz Method
    for (int t = g_data.num_time_steps - 1; t > 0; t--) {
        // Discount future cash flows back one step
        for (int i = 0; i < g_data.num_simulations; i++) {
            g_data.cash_flows[i] *= discount_factor;
        }

        // Identify in-the-money paths
        double *x_itm = (double *)malloc(g_data.num_simulations * sizeof(double));
        double *y_itm = (double *)malloc(g_data.num_simulations * sizeof(double));
        int itm_count = 0;
        for (int i = 0; i < g_data.num_simulations; i++) {
            double stock_price = g_data.asset_prices[i][t];
            if (g_data.K - stock_price > 0) {
                x_itm[itm_count] = stock_price;
                y_itm[itm_count] = g_data.cash_flows[i];
                itm_count++;
            }
        }

        if (itm_count > d) {
            // Perform least-squares regression to estimate continuation value
            int n = d + 1;
            double *beta = (double *)malloc(n * sizeof(double));
            double **A = (double **)malloc(n * sizeof(double *));
            for(int i=0; i<n; i++) A[i] = (double *)malloc(n * sizeof(double));
            double *b = (double *)malloc(n * sizeof(double));

            for (int i = 0; i < n; i++) {
                b[i] = 0.0;
                 for (int j = 0; j < n; j++) A[i][j] = 0.0;
            }

            for (int k = 0; k < itm_count; k++) {
                double x_val = x_itm[k];
                double y_val = y_itm[k];
                double p_x[n];
                p_x[0] = 1.0;
                for (int i = 1; i < n; i++) p_x[i] = p_x[i-1] * x_val;

                for (int i = 0; i < n; i++) {
                    b[i] += p_x[i] * y_val;
                    for (int j = 0; j < n; j++) {
                        A[i][j] += p_x[i] * p_x[j];
                    }
                }
            }

            if (solve_linear_system(A, b, beta, n) == 0) {
                 // Update cash flows based on optimal exercise decision
                 int current_itm_idx = 0;
                 for (int i = 0; i < g_data.num_simulations; i++) {
                     double stock_price = g_data.asset_prices[i][t];
                     if (g_data.K - stock_price > 0) {
                         double continuation_value = 0.0;
                         double x_pow = 1.0;
                         for (int j = 0; j < n; j++) {
                            continuation_value += beta[j] * x_pow;
                            x_pow *= stock_price;
                         }
                         double exercise_value = g_data.K - stock_price;
                         if (exercise_value > continuation_value) {
                             g_data.cash_flows[i] = exercise_value;
                         }
                     }
                 }
            }

            free(beta);
            for(int i=0; i<n; i++) free(A[i]);
            free(A);
            free(b);
        }
        
        free(x_itm);
        free(y_itm);
    }

    // Calculate the final option price as the discounted average of cash flows at t=1
    double total_payoff = 0.0;
    for (int i = 0; i < g_data.num_simulations; i++) {
        total_payoff += g_data.cash_flows[i];
    }
    g_data.option_price = (total_payoff / g_data.num_simulations) * discount_factor;
}

void cleanup() {
    for (int i = 0; i < g_data.num_simulations; i++) {
        free(g_data.asset_prices[i]);
    }
    free(g_data.asset_prices);
    free(g_data.cash_flows);
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
    printf("%.8f\n", g_data.option_price);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
