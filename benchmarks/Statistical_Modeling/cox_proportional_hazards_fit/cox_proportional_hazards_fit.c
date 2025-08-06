#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

// --- Mersenne Twister (MT19937) Generator --- (DO NOT MODIFY)
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

// Benchmark parameters
int NUM_SUBJECTS;
int NUM_PREDICTORS;
int NUM_ITERATIONS;

// Data arrays and model coefficients
double *X;          // Covariate matrix (NUM_SUBJECTS x NUM_PREDICTORS)
double *time_val;   // Survival times for each subject
int *status;        // Event status (1=event, 0=censored)
double *beta;       // Model coefficients

// Workspace arrays for computation
double *score;                      // Score vector (gradient)
double *risk_numerator_sum_X_buffer; // Temporary buffer

// Final result accumulator
double final_result;

// Separates data setup from the core computation.
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_subjects num_predictors num_iterations seed\n", argv[0]);
        exit(1);
    }

    NUM_SUBJECTS = atoi(argv[1]);
    NUM_PREDICTORS = atoi(argv[2]);
    NUM_ITERATIONS = atoi(argv[3]);
    uint32_t seed = atoi(argv[4]);
    mt_seed(seed);

    // Allocate memory
    X = (double *)malloc(NUM_SUBJECTS * NUM_PREDICTORS * sizeof(double));
    time_val = (double *)malloc(NUM_SUBJECTS * sizeof(double));
    status = (int *)malloc(NUM_SUBJECTS * sizeof(int));
    beta = (double *)malloc(NUM_PREDICTORS * sizeof(double));
    score = (double *)malloc(NUM_PREDICTORS * sizeof(double));
    risk_numerator_sum_X_buffer = (double *)malloc(NUM_PREDICTORS * sizeof(double));

    if (!X || !time_val || !status || !beta || !score || !risk_numerator_sum_X_buffer) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Initialize data with random values
    for (int i = 0; i < NUM_SUBJECTS; ++i) {
        time_val[i] = 1.0 + 1000.0 * ((double)mt_rand() / (double)UINT32_MAX);
        // ~80% events, ~20% censored
        status[i] = (mt_rand() % 5 != 0);
        for (int j = 0; j < NUM_PREDICTORS; ++j) {
            X[i * NUM_PREDICTORS + j] = ((double)mt_rand() / (double)UINT32_MAX) - 0.5;
        }
    }

    // Initialize model coefficients to zero
    for (int i = 0; i < NUM_PREDICTORS; ++i) {
        beta[i] = 0.0;
    }
}

// Performs the core computation of the benchmark.
void run_computation() {
    const double learning_rate = 0.01;

    for (int iter = 0; iter < NUM_ITERATIONS; ++iter) {
        // Reset score vector for this iteration
        for (int k = 0; k < NUM_PREDICTORS; ++k) {
            score[k] = 0.0;
        }

        // Loop through each subject that experienced an event
        for (int i = 0; i < NUM_SUBJECTS; ++i) {
            if (status[i] == 0) {
                continue; // Skip censored subjects
            }

            // Calculate risk set for subject i and related sums
            double risk_denominator = 0.0;
            for (int k = 0; k < NUM_PREDICTORS; ++k) {
                risk_numerator_sum_X_buffer[k] = 0.0;
            }

            // The risk set for subject i includes all subjects j with time[j] >= time[i]
            for (int j = 0; j < NUM_SUBJECTS; ++j) {
                if (time_val[j] >= time_val[i]) {
                    // Calculate linear predictor eta_j = beta * X_j
                    double eta_j = 0.0;
                    for (int k = 0; k < NUM_PREDICTORS; ++k) {
                        eta_j += beta[k] * X[j * NUM_PREDICTORS + k];
                    }
                    double exp_eta_j = exp(eta_j);

                    risk_denominator += exp_eta_j;
                    for (int k = 0; k < NUM_PREDICTORS; ++k) {
                        risk_numerator_sum_X_buffer[k] += X[j * NUM_PREDICTORS + k] * exp_eta_j;
                    }
                }
            }

            // Update score vector component for subject i
            if (risk_denominator > 1e-9) { // Avoid division by zero
                for (int k = 0; k < NUM_PREDICTORS; ++k) {
                    score[k] += X[i * NUM_PREDICTORS + k] - (risk_numerator_sum_X_buffer[k] / risk_denominator);
                }
            }
        }

        // Simple gradient ascent update for beta coefficients
        for (int j = 0; j < NUM_PREDICTORS; ++j) {
            beta[j] += learning_rate * score[j] / NUM_SUBJECTS;
        }
    }

    // Accumulate final result to prevent dead code elimination
    final_result = 0.0;
    for (int j = 0; j < NUM_PREDICTORS; ++j) {
        final_result += beta[j];
    }
}

// Frees all allocated memory.
void cleanup() {
    free(X);
    free(time_val);
    free(status);
    free(beta);
    free(score);
    free(risk_numerator_sum_X_buffer);
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
    printf("%f\n", final_result);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
