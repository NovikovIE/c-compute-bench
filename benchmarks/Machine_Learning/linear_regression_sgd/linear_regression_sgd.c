#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// Mersenne Twister (MT19937) generator - DO NOT MODIFY
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
// End of Mersenne Twister

// Benchmark parameters and data
int num_samples;
int num_features;
int num_epochs;
double learning_rate;

double *X; // Input features (num_samples x num_features)
double *y; // Target values (num_samples)
double *weights; // Model weights (num_features)
double final_result; // To prevent dead code elimination

// Function to generate a random double between 0.0 and 1.0
double rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <num_samples> <num_features> <num_epochs> <learning_rate> <seed>\n", argv[0]);
        exit(1);
    }

    num_samples = atoi(argv[1]);
    num_features = atoi(argv[2]);
    num_epochs = atoi(argv[3]);
    learning_rate = atof(argv[4]);
    uint32_t seed = (uint32_t)strtoul(argv[5], NULL, 10);

    mt_seed(seed);

    // Allocate memory
    X = (double*)malloc((size_t)num_samples * num_features * sizeof(double));
    y = (double*)malloc((size_t)num_samples * sizeof(double));
    weights = (double*)malloc((size_t)num_features * sizeof(double));

    if (!X || !y || !weights) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // To make the learning task realistic, we generate data from a known model
    // with some added noise.
    double *true_weights = (double*)malloc((size_t)num_features * sizeof(double));
    if (!true_weights) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    for (int j = 0; j < num_features; ++j) {
        true_weights[j] = rand_double() * 2.0 - 1.0; // Weights in [-1, 1]
    }

    // Generate training data
    for (int i = 0; i < num_samples; ++i) {
        double y_true = 0.0;
        for (int j = 0; j < num_features; ++j) {
            X[(size_t)i * num_features + j] = rand_double(); // Features in [0, 1]
            y_true += true_weights[j] * X[(size_t)i * num_features + j];
        }
        // Add small random noise to the target value
        y[i] = y_true + (rand_double() - 0.5) * 0.1;
    }
    
    free(true_weights); // No longer needed

    // Initialize model weights to zero
    for (int j = 0; j < num_features; ++j) {
        weights[j] = 0.0;
    }
}

void run_computation() {
    for (int epoch = 0; epoch < num_epochs; ++epoch) {
        for (int i = 0; i < num_samples; ++i) {
            // Predict y_pred = dot(weights, X_i)
            double y_pred = 0.0;
            for (int j = 0; j < num_features; ++j) {
                y_pred += weights[j] * X[(size_t)i * num_features + j];
            }

            // Calculate error
            double error = y_pred - y[i];

            // Update weights using gradient descent
            for (int j = 0; j < num_features; ++j) {
                weights[j] -= learning_rate * error * X[(size_t)i * num_features + j];
            }
        }
    }

    // Calculate a final result to prevent dead code elimination and provide a checksum.
    // The sum of the final weights is a good candidate.
    double weight_sum = 0.0;
    for (int j = 0; j < num_features; ++j) {
        weight_sum += weights[j];
    }
    final_result = weight_sum;
}

void cleanup() {
    free(X);
    free(y);
    free(weights);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%f\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
