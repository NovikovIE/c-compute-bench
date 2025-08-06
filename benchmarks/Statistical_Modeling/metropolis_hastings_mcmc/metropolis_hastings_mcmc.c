#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <string.h>

// --- MERSENNE TWISTER (DO NOT MODIFY) ---
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
int num_samples;
int num_burn_in;
int num_parameters;
int num_data_points;

// Data structures
double *data_X;          // Predictor variables: [num_data_points x num_parameters]
double *data_y;          // Response variable: [num_data_points]
double *samples;         // Stored MCMC samples: [num_samples x num_parameters]
double *current_params;  // Current state of the MCMC chain: [num_parameters]
double *proposal_params; // Proposed state for the chain: [num_parameters]

// Final result accumulator
double final_result = 0.0;

// Generates a random double in [min, max)
double random_double(double min, double max) {
    return min + ((double)mt_rand() / (double)(UINT32_MAX)) * (max - min);
}

// Calculates the negative log-posterior of a Bayesian linear regression model.
// We use a simple model: y ~ N(X*beta, 1) with prior beta ~ N(0, 1).
// The negative log-posterior is proportional to the sum of squared errors plus the L2 norm of beta.
double calculate_neg_log_posterior(const double* params) {
    double log_likelihood = 0.0;

    // Likelihood part: Sum of Squared Errors
    for (int i = 0; i < num_data_points; ++i) {
        double prediction = 0.0;
        for (int j = 0; j < num_parameters; ++j) {
            prediction += data_X[i * num_parameters + j] * params[j];
        }
        double error = data_y[i] - prediction;
        log_likelihood += error * error; // Simplified, assuming variance=1
    }

    // Prior part: L2 regularization on parameters
    double log_prior = 0.0;
    for (int j = 0; j < num_parameters; ++j) {
        log_prior += params[j] * params[j]; // Simplified, assuming variance=1
    }

    return 0.5 * (log_likelihood + log_prior);
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s num_samples num_burn_in num_parameters num_data_points seed\n", argv[0]);
        exit(1);
    }

    num_samples = atoi(argv[1]);
    num_burn_in = atoi(argv[2]);
    num_parameters = atoi(argv[3]);
    num_data_points = atoi(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    mt_seed(seed);

    // Allocate memory
    data_X = (double*)malloc(num_data_points * num_parameters * sizeof(double));
    data_y = (double*)malloc(num_data_points * sizeof(double));
    samples = (double*)malloc(num_samples * num_parameters * sizeof(double));
    current_params = (double*)malloc(num_parameters * sizeof(double));
    proposal_params = (double*)malloc(num_parameters * sizeof(double));
    double *true_beta = (double*)malloc(num_parameters * sizeof(double));

    if (!data_X || !data_y || !samples || !current_params || !proposal_params || !true_beta) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(1);
    }

    // Generate synthetic data
    // 1. Generate "true" parameters
    for (int j = 0; j < num_parameters; ++j) {
        true_beta[j] = random_double(-2.0, 2.0);
    }

    // 2. Generate predictor data X and response data y
    for (int i = 0; i < num_data_points; ++i) {
        data_y[i] = random_double(-0.5, 0.5); // Add noise term first
        for (int j = 0; j < num_parameters; ++j) {
            data_X[i * num_parameters + j] = random_double(-1.0, 1.0);
            data_y[i] += data_X[i * num_parameters + j] * true_beta[j];
        }
    }

    // 3. Initialize the MCMC chain
    for (int j = 0; j < num_parameters; ++j) {
        current_params[j] = random_double(-0.1, 0.1);
    }

    free(true_beta); // No longer needed
}

void run_computation() {
    int total_iterations = num_samples + num_burn_in;
    double proposal_width = 0.05;
    int sample_count = 0;

    double neg_log_post_current = calculate_neg_log_posterior(current_params);

    for (int i = 0; i < total_iterations; ++i) {
        // 1. Propose new parameters (symmetric proposal)
        for (int j = 0; j < num_parameters; ++j) {
            proposal_params[j] = current_params[j] + random_double(-proposal_width, proposal_width);
        }

        // 2. Calculate acceptance probability in log space
        double neg_log_post_proposal = calculate_neg_log_posterior(proposal_params);
        double log_acceptance_ratio = neg_log_post_current - neg_log_post_proposal;

        // 3. Accept or reject the proposal
        if (log(random_double(0.0, 1.0)) < log_acceptance_ratio) {
            memcpy(current_params, proposal_params, num_parameters * sizeof(double));
            neg_log_post_current = neg_log_post_proposal;
        }

        // 4. Store the sample if past the burn-in period
        if (i >= num_burn_in) {
            memcpy(&samples[sample_count * num_parameters], current_params, num_parameters * sizeof(double));
            sample_count++;
        }
    }

    // Accumulate a result to prevent dead code elimination
    for (int i = 0; i < num_samples; ++i) {
        for (int j = 0; j < num_parameters; ++j) {
            final_result += samples[i * num_parameters + j];
        }
    }
}

void cleanup() {
    free(data_X);
    free(data_y);
    free(samples);
    free(current_params);
    free(proposal_params);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%f\n", final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
