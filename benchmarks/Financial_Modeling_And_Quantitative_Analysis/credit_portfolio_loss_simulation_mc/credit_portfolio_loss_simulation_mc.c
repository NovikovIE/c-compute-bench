#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- Mersenne Twister (MT19937) Generator (Verbatim) ---
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
// --- End of MT19937 ---

// --- Benchmark Globals ---
static int num_obligors;
static int num_simulations;
static double correlation_factor;
static double default_threshold;

// Pre-computed random factors to isolate simulation computation
static double* systemic_factors;      // Z for each simulation
static double* idiosyncratic_factors; // E_i for each obligor in each simulation

// Final result
static unsigned long long total_defaults = 0;

#define PI 3.14159265358979323846

// --- Utility Function ---
/**
 * @brief Generates a standard normal random variate using the Box-Muller transform.
 * 
 * This function uses two uniform random numbers to generate two normal random numbers.
 * A static variable is used to cache the second value for the next call, improving efficiency.
 * @return A double drawn from a standard normal distribution N(0,1).
 */
static double generate_normal() {
    static int has_cached_value = 0;
    static double cached_value;

    if (has_cached_value) {
        has_cached_value = 0;
        return cached_value;
    }

    double u1, u2;
    do {
        u1 = (double)mt_rand() / UINT32_MAX;
    } while (u1 == 0.0); // Avoid log(0)
    u2 = (double)mt_rand() / UINT32_MAX;

    double r = sqrt(-2.0 * log(u1));
    double theta = 2.0 * PI * u2;

    cached_value = r * sin(theta);
    has_cached_value = 1;

    return r * cos(theta);
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <num_obligors> <num_simulations> <correlation_factor> <default_threshold> <seed>\n", argv[0]);
        exit(1);
    }

    num_obligors = atoi(argv[1]);
    num_simulations = atoi(argv[2]);
    correlation_factor = atof(argv[3]);
    default_threshold = atof(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    mt_seed(seed);

    systemic_factors = (double *)malloc(num_simulations * sizeof(double));
    if (!systemic_factors) {
        perror("Failed to allocate memory for systemic_factors");
        exit(1);
    }

    size_t num_idiosyncratic_factors = (size_t)num_simulations * num_obligors;
    idiosyncratic_factors = (double *)malloc(num_idiosyncratic_factors * sizeof(double));
    if (!idiosyncratic_factors) {
        perror("Failed to allocate memory for idiosyncratic_factors");
        free(systemic_factors);
        exit(1);
    }

    for (int i = 0; i < num_simulations; ++i) {
        systemic_factors[i] = generate_normal();
    }

    for (size_t i = 0; i < num_idiosyncratic_factors; ++i) {
        idiosyncratic_factors[i] = generate_normal();
    }
}

void run_computation() {
    unsigned long long cumulative_defaults = 0;
    const double sqrt_corr = sqrt(correlation_factor);
    const double sqrt_one_minus_corr = sqrt(1.0 - correlation_factor);
    const size_t num_obligors_size_t = num_obligors;

    for (int s = 0; s < num_simulations; ++s) {
        const double z = systemic_factors[s];
        const size_t offset = (size_t)s * num_obligors_size_t;
        for (size_t o = 0; o < num_obligors_size_t; ++o) {
            const double e = idiosyncratic_factors[offset + o];
            const double asset_value = sqrt_corr * z + sqrt_one_minus_corr * e;
            if (asset_value < default_threshold) {
                cumulative_defaults++;
            }
        }
    }
    total_defaults = cumulative_defaults;
}

void cleanup() {
    free(systemic_factors);
    free(idiosyncratic_factors);
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
    printf("%llu\n", total_defaults);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
