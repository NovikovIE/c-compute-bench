#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>
#include <string.h>

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

// Benchmark parameters
static int num_data_points;
static int num_parameters;
static int max_iterations;
static double learning_rate;
static uint32_t seed;

// Data structures
static double* data_X;       // Feature matrix (flattened)
static int* data_y;          // Target labels
static double* parameters_w; // Model parameters (weights)
static double* gradient;     // Gradient vector

// Final result
static double final_result;

// Helper to generate a random double between 0.0 and 1.0
double rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <num_data_points> <num_parameters> <max_iterations> <learning_rate> <seed>\n", argv[0]);
        exit(1);
    }

    num_data_points = atoi(argv[1]);
    num_parameters = atoi(argv[2]);
    max_iterations = atoi(argv[3]);
    learning_rate = atof(argv[4]);
    seed = (uint32_t)atoi(argv[5]);

    mt_seed(seed);

    // Allocate memory
    data_X = (double*)malloc(num_data_points * num_parameters * sizeof(double));
    data_y = (int*)malloc(num_data_points * sizeof(int));
    parameters_w = (double*)malloc(num_parameters * sizeof(double));
    gradient = (double*)malloc(num_parameters * sizeof(double));

    if (!data_X || !data_y || !parameters_w || !gradient) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize data with random values
    for (int i = 0; i < num_data_points; ++i) {
        for (int j = 0; j < num_parameters; ++j) {
            data_X[i * num_parameters + j] = rand_double() - 0.5;
        }
        data_y[i] = mt_rand() % 2; // Binary classification labels (0 or 1)
    }

    // Initialize parameters to zero
    memset(parameters_w, 0, num_parameters * sizeof(double));
}

void run_computation() {
    // Gradient Ascent for Logistic Regression (MLE)
    for (int iter = 0; iter < max_iterations; ++iter) {
        // Reset gradient for this iteration
        memset(gradient, 0, num_parameters * sizeof(double));

        // Calculate gradient over all data points
        for (int i = 0; i < num_data_points; ++i) {
            // Calculate dot product: w^T * x
            double dot_product = 0.0;
            for (int j = 0; j < num_parameters; ++j) {
                dot_product += parameters_w[j] * data_X[i * num_parameters + j];
            }

            // Apply sigmoid function to get probability
            // Clamp dot_product to avoid overflow in exp()
            if (dot_product > 30.0) dot_product = 30.0;
            if (dot_product < -30.0) dot_product = -30.0;
            double probability = 1.0 / (1.0 + exp(-dot_product));

            // Calculate error (y - p)
            double error = (double)data_y[i] - probability;

            // Accumulate gradient: gradient += error * x
            for (int j = 0; j < num_parameters; ++j) {
                gradient[j] += error * data_X[i * num_parameters + j];
            }
        }

        // Update parameters: w = w + learning_rate * (gradient / N)
        for (int j = 0; j < num_parameters; ++j) {
            parameters_w[j] += learning_rate * gradient[j] / (double)num_data_points;
        }
    }

    // Calculate a final result to prevent dead code elimination
    final_result = 0.0;
    for (int j = 0; j < num_parameters; ++j) {
        final_result += parameters_w[j];
    }
}

void cleanup() {
    free(data_X);
    free(data_y);
    free(parameters_w);
    free(gradient);
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
