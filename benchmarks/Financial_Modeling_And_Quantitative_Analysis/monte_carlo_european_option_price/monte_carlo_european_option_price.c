#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// For compatibility, define M_PI if not defined
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- BEGIN MERSENNE TWISTER (DO NOT MODIFY) ---
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
// --- END MERSENNE TWISTER ---

// Benchmark parameters
int num_simulations;
int num_time_steps;
int num_assets;

// Financial model parameters
double* s0;     // Initial asset prices
double* k;      // Strike prices
double* r;      // Risk-free interest rates
double* sigma;  // Volatilities
double t_maturity; // Time to maturity (in years)

// Output data
double* option_prices;
double final_result; // Accumulated result to prevent dead-code elimination

// Utility function to generate a random double between a and b
double rand_double(double a, double b) {
    return a + (mt_rand() / (double)UINT32_MAX) * (b - a);
}

// Utility function to generate a standard normal random variable using Box-Muller transform
double generate_normal() {
    // Ensure u1 is not 0 to avoid log(0)
    double u1 = (mt_rand() + 1.0) / (UINT32_MAX + 2.0);
    double u2 = mt_rand() / (double)UINT32_MAX;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_simulations> <num_time_steps> <num_assets> <seed>\n", argv[0]);
        exit(1);
    }

    num_simulations = atoi(argv[1]);
    num_time_steps = atoi(argv[2]);
    num_assets = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    mt_seed(seed);

    // Allocate memory on the heap
    s0 = (double*)malloc(num_assets * sizeof(double));
    k = (double*)malloc(num_assets * sizeof(double));
    r = (double*)malloc(num_assets * sizeof(double));
    sigma = (double*)malloc(num_assets * sizeof(double));
    option_prices = (double*)malloc(num_assets * sizeof(double));

    if (!s0 || !k || !r || !sigma || !option_prices) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Common time to maturity for all options
    t_maturity = 1.0; 

    // Generate random financial parameters for each asset
    for (int i = 0; i < num_assets; i++) {
        s0[i] = rand_double(90.0, 110.0);   // Initial price around 100
        k[i] = rand_double(95.0, 105.0);    // Strike price near the initial price
        r[i] = rand_double(0.01, 0.05);     // Risk-free rate between 1% and 5%
        sigma[i] = rand_double(0.1, 0.4);   // Volatility between 10% and 40%
    }
    
    final_result = 0.0;
}

void run_computation() {
    double dt = t_maturity / (double)num_time_steps;
    double sqrt_dt = sqrt(dt);

    for (int i = 0; i < num_assets; i++) {
        double total_payoff = 0.0;
        double current_s0 = s0[i];
        double current_k = k[i];
        double current_r = r[i];
        double current_sigma = sigma[i];
        double drift = (current_r - 0.5 * current_sigma * current_sigma) * dt;
        double vol_term = current_sigma * sqrt_dt;

        for (int j = 0; j < num_simulations; j++) {
            double asset_price = current_s0;
            for (int step = 0; step < num_time_steps; step++) {
                double z = generate_normal();
                asset_price *= exp(drift + vol_term * z);
            }

            // Calculate payoff for a European call option
            double payoff = asset_price > current_k ? asset_price - current_k : 0.0;
            total_payoff += payoff;
        }

        double avg_payoff = total_payoff / (double)num_simulations;
        // Discount the average payoff back to present value
        option_prices[i] = exp(-current_r * t_maturity) * avg_payoff;
    }

    // Accumulate results to prevent dead-code elimination
    double sum = 0.0;
    for (int i = 0; i < num_assets; i++) {
        sum += option_prices[i];
    }
    final_result = sum;
}

void cleanup() {
    free(s0);
    free(k);
    free(r);
    free(sigma);
    free(option_prices);
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
    printf("%f\n", final_result);

    // Print the timing info to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
