#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

// --- BEGIN MERSENNE TWISTER (MT19937) ---
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
typedef struct {
    int time_series_length;
    int p_order; // ARCH component order
    int q_order; // GARCH component order
    double* returns;
    double* alpha;   // ARCH parameters (size p)
    double* beta;    // GARCH parameters (size q)
    double* sigma2;  // Conditional variance series
    double omega;   // Constant term
    double final_log_likelihood;
} GarchData;

GarchData g_data;

// Helper to generate a random double in [0.0, 1.0]
double mt_rand_double(void) {
    return (double)mt_rand() / 4294967295.0; // 2^32 - 1
}

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <time_series_length> <p_order> <q_order> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.time_series_length = atoi(argv[1]);
    g_data.p_order = atoi(argv[2]);
    g_data.q_order = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    mt_seed(seed);

    if (g_data.time_series_length <= 0 || g_data.p_order <= 0 || g_data.q_order <= 0) {
        fprintf(stderr, "FATAL: All parameters must be positive integers.\n");
        exit(1);
    }
    if (g_data.p_order >= g_data.time_series_length || g_data.q_order >= g_data.time_series_length) {
        fprintf(stderr, "FATAL: p_order and q_order must be smaller than time_series_length.\n");
        exit(1);
    }

    // Allocate memory on the heap
    g_data.returns = (double*)malloc(g_data.time_series_length * sizeof(double));
    g_data.alpha = (double*)malloc(g_data.p_order * sizeof(double));
    g_data.beta = (double*)malloc(g_data.q_order * sizeof(double));
    g_data.sigma2 = (double*)malloc(g_data.time_series_length * sizeof(double));

    if (!g_data.returns || !g_data.alpha || !g_data.beta || !g_data.sigma2) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Generate synthetic financial returns data (e.g., daily returns)
    for (int i = 0; i < g_data.time_series_length; ++i) {
        g_data.returns[i] = mt_rand_double() * 0.1 - 0.05; // Returns between -5% and +5%
    }

    // Initialize GARCH model parameters with small positive random values
    g_data.omega = mt_rand_double() * 0.00001;
    double alpha_sum = 0.0;
    for (int i = 0; i < g_data.p_order; ++i) {
        g_data.alpha[i] = mt_rand_double() * 0.1; 
        alpha_sum += g_data.alpha[i];
    }
    double beta_sum = 0.0;
    for (int i = 0; i < g_data.q_order; ++i) {
        g_data.beta[i] = mt_rand_double() * 0.8; 
        beta_sum += g_data.beta[i];
    }

    // Enforce stationarity condition (alpha + beta < 1) by scaling if necessary
    if (alpha_sum + beta_sum >= 1.0) {
        double scale = 0.99 / (alpha_sum + beta_sum);
        for (int i = 0; i < g_data.p_order; ++i) g_data.alpha[i] *= scale;
        for (int i = 0; i < g_data.q_order; ++i) g_data.beta[i] *= scale;
    }
}

void run_computation() {
    const int NUM_ITERATIONS = 25; // Simulate steps of a numerical optimization
    const double LOG2PI = log(2.0 * 3.14159265358979323846);

    int p = g_data.p_order;
    int q = g_data.q_order;
    int N = g_data.time_series_length;
    int start_offset = (p > q) ? p : q;

    double total_log_likelihood = 0.0;
    
    // Calculate initial variance once using the whole series for a stable starting point
    double initial_var = 0.0;
    for(int i = 0; i < N; ++i) {
        initial_var += g_data.returns[i] * g_data.returns[i];
    }
    initial_var /= N;
    if (initial_var < 1e-9) initial_var = 1e-9; // Avoid zero variance

    for (int iter = 0; iter < NUM_ITERATIONS; ++iter) {
        total_log_likelihood = 0.0;

        // Initialize early sigma2 values for the GARCH recursion
        for (int i = 0; i < start_offset; ++i) {
            g_data.sigma2[i] = initial_var;
        }

        // Main loop to calculate conditional variance and log-likelihood
        for (int t = start_offset; t < N; ++t) {
            double arch_sum = 0.0;
            for (int i = 1; i <= p; ++i) {
                arch_sum += g_data.alpha[i - 1] * g_data.returns[t - i] * g_data.returns[t - i];
            }

            double garch_sum = 0.0;
            for (int j = 1; j <= q; ++j) {
                garch_sum += g_data.beta[j - 1] * g_data.sigma2[t - j];
            }

            g_data.sigma2[t] = g_data.omega + arch_sum + garch_sum;
            
            if (g_data.sigma2[t] < 1e-9) g_data.sigma2[t] = 1e-9; // Ensure positivity
            
            // Log-likelihood for a normal distribution: -0.5 * (log(2*pi) + log(sigma^2) + e^2/sigma^2)
            total_log_likelihood -= 0.5 * (LOG2PI + log(g_data.sigma2[t]) + (g_data.returns[t] * g_data.returns[t]) / g_data.sigma2[t]);
        }
        
        // Perturb a parameter slightly to ensure the loop is not optimized away by the compiler
        g_data.omega *= 1.00000001;
    }

    g_data.final_log_likelihood = total_log_likelihood;
}

void cleanup() {
    free(g_data.returns);
    free(g_data.alpha);
    free(g_data.beta);
    free(g_data.sigma2);
}

// --- MAIN FUNCTION ---
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    cleanup();

    // Print the final computational result to stdout
    printf("%f\n", g_data.final_log_likelihood);

    // Print the timing information to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
