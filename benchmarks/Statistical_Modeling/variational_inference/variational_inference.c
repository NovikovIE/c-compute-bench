#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// Use a constant for PI to be portable
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- START MERSENNE TWISTER (MT19937) ---
// This implementation is provided and should not be modified.
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
// --- END MERSENNE TWISTER ---

// Helper to generate a random number from a standard normal distribution
// using the Box-Muller transform.
static double rand_normal(double mean, double stddev) {
    static int have_spare = 0;
    static double spare;

    if (have_spare) {
        have_spare = 0;
        return mean + stddev * spare;
    }

    have_spare = 1;
    double u1, u2;
    // Ensure u1 is not zero to avoid log(0)
    do {
        u1 = (double)mt_rand() / UINT32_MAX;
    } while (u1 == 0.0);
    u2 = (double)mt_rand() / UINT32_MAX;
    
    double mag = sqrt(-2.0 * log(u1));
    spare = mag * sin(2.0 * M_PI * u2);
    return mean + stddev * (mag * cos(2.0 * M_PI * u2));
}

// Benchmark parameters
static int num_data_points;
static int num_parameters; // Represents the number of Gaussian components (K)
static int num_iterations;

// Data and model parameters
static double *data_points;      // [num_data_points]
static double *means;            // [num_parameters]
static double *variances;        // [num_parameters]
static double *weights;          // [num_parameters]
static double *responsibilities; // [num_data_points * num_parameters]

// Final result to prevent dead code elimination
static double final_result = 0.0;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_data_points num_parameters num_iterations seed\n", argv[0]);
        exit(1);
    }

    num_data_points = strtol(argv[1], NULL, 10);
    num_parameters = strtol(argv[2], NULL, 10);
    num_iterations = strtol(argv[3], NULL, 10);
    uint32_t seed = (uint32_t)strtoul(argv[4], NULL, 10);

    if (num_data_points <= 0 || num_parameters <= 0 || num_iterations <= 0) {
        fprintf(stderr, "Parameters must be positive integers.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate memory
    data_points = (double *)malloc(num_data_points * sizeof(double));
    means = (double *)malloc(num_parameters * sizeof(double));
    variances = (double *)malloc(num_parameters * sizeof(double));
    weights = (double *)malloc(num_parameters * sizeof(double));
    responsibilities = (double *)malloc(num_data_points * num_parameters * sizeof(double));

    if (!data_points || !means || !variances || !weights || !responsibilities) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Generate synthetic data from a known mixture of Gaussians
    for (int i = 0; i < num_data_points; ++i) {
        int component = mt_rand() % num_parameters;
        double mean = component * 50.0;
        double stddev = 5.0 + (component * 2.0);
        data_points[i] = rand_normal(mean, stddev);
    }

    // Initialize model parameters
    // - Weights: Uniform
    // - Means: Randomly chosen data points
    // - Variances: Global variance of the data (or simply 1.0 for simplicity)
    double total_variance = 0;
    double total_mean = 0;
    for(int i = 0; i < num_data_points; ++i) total_mean += data_points[i];
    total_mean /= num_data_points;
    for(int i = 0; i < num_data_points; ++i) total_variance += (data_points[i]-total_mean)*(data_points[i]-total_mean);
    total_variance /= num_data_points;

    for (int k = 0; k < num_parameters; ++k) {
        weights[k] = 1.0 / num_parameters;
        means[k] = data_points[mt_rand() % num_data_points];
        variances[k] = total_variance;
    }
}

void run_computation() {
    double current_log_likelihood = 0.0;
    
    for (int iter = 0; iter < num_iterations; ++iter) {
        // --- E-step: Calculate responsibilities --- 
        current_log_likelihood = 0.0;
        for (int n = 0; n < num_data_points; ++n) {
            double denom = 0.0;
            for (int k = 0; k < num_parameters; ++k) {
                double var_k = variances[k];
                if (var_k <= 0) var_k = 1e-9; // safety check
                double exponent = -0.5 * (data_points[n] - means[k]) * (data_points[n] - means[k]) / var_k;
                double pdf = (1.0 / sqrt(2.0 * M_PI * var_k)) * exp(exponent);
                double weighted_pdf = weights[k] * pdf;
                responsibilities[n * num_parameters + k] = weighted_pdf;
                denom += weighted_pdf;
            }
            
            current_log_likelihood += log(denom > 1e-9 ? denom : 1e-9);

            for (int k = 0; k < num_parameters; ++k) {
                 responsibilities[n * num_parameters + k] /= (denom > 1e-9 ? denom : 1.0);
            }
        }

        // --- M-step: Update model parameters --- 
        for (int k = 0; k < num_parameters; ++k) {
            double sum_resp = 0.0;
            for (int n = 0; n < num_data_points; ++n) {
                sum_resp += responsibilities[n * num_parameters + k];
            }

            // Skip update if component has collapsed
            if (sum_resp < 1e-9) {
                continue;
            }

            // Update weight
            weights[k] = sum_resp / num_data_points;
            
            // Update mean
            double new_mean = 0.0;
            for (int n = 0; n < num_data_points; ++n) {
                new_mean += responsibilities[n * num_parameters + k] * data_points[n];
            }
            means[k] = new_mean / sum_resp;

            // Update variance
            double new_var = 0.0;
            for (int n = 0; n < num_data_points; ++n) {
                new_var += responsibilities[n * num_parameters + k] * (data_points[n] - means[k]) * (data_points[n] - means[k]);
            }
            variances[k] = (new_var / sum_resp) + 1e-6; // Add epsilon for stability
        }
    }

    final_result = current_log_likelihood;
}

void cleanup() {
    free(data_points);
    free(means);
    free(variances);
    free(weights);
    free(responsibilities);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final log-likelihood to stdout
    printf("%f\n", final_result);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
