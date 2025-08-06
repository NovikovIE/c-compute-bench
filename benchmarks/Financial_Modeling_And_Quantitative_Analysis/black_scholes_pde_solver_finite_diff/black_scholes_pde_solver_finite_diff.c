#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define MAX(a, b) (((a) > (b)) ? (a) : (b))

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
// --- END: Mersenne Twister ---

// Benchmark parameters
int GRID_SIZE_PRICE; // Number of price steps (N)
int GRID_SIZE_TIME;  // Number of time steps (M)

// Financial model parameters
double S_MAX; // Max asset price
double K;     // Strike price
double R;     // Risk-free rate
double SIGMA; // Volatility
double T;     // Time to maturity (in years)

// Data grids for computation
double* grid_a;
double* grid_b;

// Final result
double final_result;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s grid_size_price grid_size_time seed\n", argv[0]);
        exit(1);
    }

    GRID_SIZE_PRICE = atoi(argv[1]);
    GRID_SIZE_TIME = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    mt_seed(seed); // Seed the generator, though not used in this deterministic model

    // Set financial parameters for a European Call Option
    S_MAX = 200.0;  // Maximum price of the underlying asset to consider
    K = 100.0;      // Strike price
    R = 0.05;       // Risk-free interest rate (5%)
    SIGMA = 0.20;   // Volatility (20%)
    T = 1.0;        // Time to maturity (1 year)

    // Allocate memory for the two price grids (current and previous time steps)
    grid_a = (double*)malloc((GRID_SIZE_PRICE + 1) * sizeof(double));
    grid_b = (double*)malloc((GRID_SIZE_PRICE + 1) * sizeof(double));

    if (!grid_a || !grid_b) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Initialize the grid at maturity (t=T)
    // For a call option, the value is max(S - K, 0)
    double dS = S_MAX / GRID_SIZE_PRICE;
    for (int i = 0; i <= GRID_SIZE_PRICE; ++i) {
        double current_price = i * dS;
        grid_a[i] = MAX(current_price - K, 0.0);
    } 
}

void run_computation() {
    double dt = T / GRID_SIZE_TIME;
    double dS = S_MAX / GRID_SIZE_PRICE;

    double* current_prices = grid_a;
    double* next_prices = grid_b;

    // Solve the PDE by iterating backwards in time from T to 0
    for (int j = GRID_SIZE_TIME - 1; j >= 0; --j) {
        // Set boundary conditions for this time step
        // At S=0, option value is 0
        next_prices[0] = 0.0;
        // At S=S_max, option value is S_max - K*exp(-r*(T-t))
        double current_time = j * dt;
        next_prices[GRID_SIZE_PRICE] = S_MAX - K * exp(-R * (T - current_time));

        // Calculate interior grid points using the explicit finite difference method
        for (int i = 1; i < GRID_SIZE_PRICE; ++i) {
            double i_double = (double)i;
            double sigma_sq_i_sq = SIGMA * SIGMA * i_double * i_double;
            double r_i = R * i_double;

            double a = 0.5 * dt * (sigma_sq_i_sq - r_i);
            double b = 1.0 - dt * (sigma_sq_i_sq + R);
            double c = 0.5 * dt * (sigma_sq_i_sq + r_i);

            next_prices[i] = a * current_prices[i - 1] + 
                               b * current_prices[i] + 
                               c * current_prices[i + 1];
        }

        // Swap pointers for the next iteration
        double* temp = current_prices;
        current_prices = next_prices;
        next_prices = temp;
    }

    // Accumulate the final results from the grid at t=0 to prevent dead code elimination
    // The result is stored in the `current_prices` array after the loops finish.
    final_result = 0.0;
    for (int i = 0; i <= GRID_SIZE_PRICE; i++) {
        final_result += current_prices[i];
    }
}

void cleanup() {
    free(grid_a);
    free(grid_b);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final accumulated result to stdout
    printf("%f\n", final_result);

    // Print timing information to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
