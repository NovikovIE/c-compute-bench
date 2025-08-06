#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

// --- MERSENNE TWISTER (Verbatim as per requirements) ---
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

// --- BENCHMARK DATA AND PARAMETERS ---

// Fixed internal parameters
#define NUM_PATHS 400

// Global struct to hold all benchmark data
typedef struct {
    int num_market_instruments;
    int time_series_length;

    // Model parameters
    double a; // Mean reversion speed
    double sigma; // Volatility
    double theta; // Long-term mean level of the rate
    double initial_short_rate; // r(0)

    // Market Data (Inputs)
    double* market_bond_prices;
    double* bond_maturities; // T for each bond

    // Calculation result
    double total_squared_error;
} BenchmarkData;

static BenchmarkData g_data;

// --- RANDOM NUMBER HELPERS ---

// Generates a random double in (0, 1)
double rand_double() {
    return ((double)mt_rand() + 1.0) / ((double)UINT32_MAX + 2.0);
}

// Generates a standard normal random variable using the Box-Muller transform
double rand_normal() {
    double u1 = rand_double();
    double u2 = rand_double();
    return sqrt(-2.0 * log(u1)) * cos(2.0 * 3.14159265358979323846 * u2);
}


// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_market_instruments> <time_series_length> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_market_instruments = atoi(argv[1]);
    g_data.time_series_length = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    mt_seed(seed);

    // Initialize Hull-White model parameters
    g_data.a = 0.1;      // Mean reversion speed
    g_data.sigma = 0.02;   // Volatility of the short rate
    g_data.theta = 0.04;   // Long-term mean level
    g_data.initial_short_rate = 0.03; // r(0)

    // Allocate memory on the heap
    g_data.market_bond_prices = (double*)malloc(g_data.num_market_instruments * sizeof(double));
    g_data.bond_maturities = (double*)malloc(g_data.num_market_instruments * sizeof(double));

    if (!g_data.market_bond_prices || !g_data.bond_maturities) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Generate synthetic market data for zero-coupon bonds
    for (int i = 0; i < g_data.num_market_instruments; ++i) {
        // Maturities from 0.25 years to N * 0.25 years
        double T = (double)(i + 1) * 0.25;
        g_data.bond_maturities[i] = T;

        // Generate a plausible 'market' price using a simplified discount factor
        // with some noise to simulate market imperfections.
        double noise = 1.0 + (rand_double() - 0.5) * 0.02; // +/- 1% noise
        g_data.market_bond_prices[i] = exp(-g_data.initial_short_rate * T) * noise;
    }

    g_data.total_squared_error = 0.0;
}

void run_computation() {
    double total_sq_err = 0.0;

    // For each market instrument, price it using the model and compare to market price.
    for (int i = 0; i < g_data.num_market_instruments; ++i) {
        double T = g_data.bond_maturities[i];
        double dt = T / (double)g_data.time_series_length;
        double sqrt_dt = sqrt(dt);

        double sum_of_path_prices = 0.0;

        // Monte Carlo simulation loop
        for (int p = 0; p < NUM_PATHS; ++p) {
            double r = g_data.initial_short_rate;
            double integral_r = 0.0;

            // Time-stepping loop to simulate one interest rate path (Euler-Maruyama)
            for (int s = 0; s < g_data.time_series_length; ++s) {
                // dW is the Wiener process increment
                double dW = rand_normal() * sqrt_dt;

                // Hull-White model: dr(t) = (theta(t) - a*r(t))dt + sigma*dW(t)
                // We use a constant theta for simplicity.
                double dr = (g_data.theta - g_data.a * r) * dt + g_data.sigma * dW;
                r += dr;

                // Integrate r(t) over the path to calculate the bond price
                integral_r += r * dt;
            }
            
            // Price of a zero-coupon bond is E[exp(-integral(r(s)ds))]
            double path_price = exp(-integral_r);
            sum_of_path_prices += path_price;
        }

        double model_price = sum_of_path_prices / (double)NUM_PATHS;
        double error = model_price - g_data.market_bond_prices[i];
        total_sq_err += error * error;
    }

    g_data.total_squared_error = total_sq_err;
}

void cleanup() {
    free(g_data.market_bond_prices);
    free(g_data.bond_maturities);
}


int main(int argc, char* argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout
    printf("%.8f\n", g_data.total_squared_error);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
