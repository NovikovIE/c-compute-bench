#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- Mersenne Twister (MT19937) Generator (verbatim) ---
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

// --- Benchmark Globals ---

// Parameters
int num_samples;
int num_burn_in;
int num_variables;

// Data arrays
double* variables;      // Current state of all variables in the model
double* samples_sum;    // Sum of samples after burn-in period

// Final result accumulator
double final_result;

// Helper to generate a normally distributed random number (Box-Muller transform)
double generate_normal_sample(double mean, double std_dev) {
    static int has_cached_sample = 0;
    static double cached_sample; // This will store a standard normal variate N(0,1)
    const double PI = 3.14159265358979323846;

    if (has_cached_sample) {
        has_cached_sample = 0;
        return mean + cached_sample * std_dev;
    }

    double u1, u2;
    do {
        u1 = (double)mt_rand() / (UINT32_MAX + 1.0);
    } while (u1 == 0.0); // Avoid log(0)
    u2 = (double)mt_rand() / (UINT32_MAX + 1.0);

    double mag = sqrt(-2.0 * log(u1));
    double z1 = mag * cos(2.0 * PI * u2); // N(0,1) sample
    double z2 = mag * sin(2.0 * PI * u2); // N(0,1) sample

    cached_sample = z2;
    has_cached_sample = 1;

    return mean + z1 * std_dev;
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_samples num_burn_in num_variables seed\n", argv[0]);
        exit(1);
    }

    num_samples = atoi(argv[1]);
    num_burn_in = atoi(argv[2]);
    num_variables = atoi(argv[3]);
    uint32_t seed = (uint32_t)strtoul(argv[4], NULL, 10);

    if (num_samples <= 0 || num_burn_in < 0 || num_variables <= 1) {
        fprintf(stderr, "Error: Invalid parameters. num_variables must be > 1.\n");
        exit(1);
    }

    mt_seed(seed);

    variables = (double*)malloc(num_variables * sizeof(double));
    samples_sum = (double*)malloc(num_variables * sizeof(double));

    if (!variables || !samples_sum) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    for (int i = 0; i < num_variables; ++i) {
        variables[i] = (double)mt_rand() / (UINT32_MAX + 1.0); // Initial state [0, 1)
        samples_sum[i] = 0.0;
    }
}

void run_computation() {
    int total_iterations = num_samples + num_burn_in;
    const double std_dev = 1.0; // Assume constant variance for simplicity

    for (int i = 0; i < total_iterations; ++i) {
        // Perform one full Gibbs scan (update each variable once)
        for (int j = 0; j < num_variables; ++j) {
            int prev_idx = (j == 0) ? num_variables - 1 : j - 1;
            int next_idx = (j == num_variables - 1) ? 0 : j + 1;

            // The conditional distribution of X_j depends on its neighbors X_{j-1} and X_{j+1}
            // We model this as N((X_{j-1} + X_{j+1}) / 2, 1.0)
            double conditional_mean = (variables[prev_idx] + variables[next_idx]) / 2.0;

            variables[j] = generate_normal_sample(conditional_mean, std_dev);
        }

        // If we are past the burn-in period, accumulate the samples
        if (i >= num_burn_in) {
            for (int j = 0; j < num_variables; ++j) {
                samples_sum[j] += variables[j];
            }
        }
    }

    // Calculate a final result to prevent dead code elimination.
    // We use the sum of all sample means.
    double total_sum_of_means = 0.0;
    for (int j = 0; j < num_variables; ++j) {
        total_sum_of_means += samples_sum[j];
    }

    final_result = (num_samples > 0) ? total_sum_of_means / num_samples : 0.0;
}

void cleanup() {
    free(variables);
    free(samples_sum);
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

    // Print the final result to stdout
    printf("%f\n", final_result);

    // Print the execution time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
