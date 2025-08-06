#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- Mersenne Twister (MT19937) Generator ---
// Do Not Modify
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
static int num_samples;
static int num_features;
static int num_iterations;

// Data structures
static float* data_points; // Flattened 2D array: num_samples x num_features
static int* labels;        // Array of labels: num_samples
static float* weights;     // SVM model weights: num_features

// Final result accumulator
static float final_result;

// Generates a random float between -1.0 and 1.0
static float rand_float(void) {
    return ((float)mt_rand() / (float)UINT32_MAX) * 2.0f - 1.0f;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_samples> <num_features> <num_iterations> <seed>\n", argv[0]);
        exit(1);
    }

    num_samples = atoi(argv[1]);
    num_features = atoi(argv[2]);
    num_iterations = atoi(argv[3]);
    uint32_t seed = (uint32_t)strtoul(argv[4], NULL, 10);

    mt_seed(seed);

    // Allocate memory on the heap
    data_points = (float*)malloc((size_t)num_samples * num_features * sizeof(float));
    labels = (int*)malloc((size_t)num_samples * sizeof(int));
    weights = (float*)malloc((size_t)num_features * sizeof(float));

    if (!data_points || !labels || !weights) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Initialize data with random values
    for (int i = 0; i < num_samples; ++i) {
        for (int j = 0; j < num_features; ++j) {
            data_points[i * num_features + j] = rand_float();
        }
        // Label is either -1 or 1
        labels[i] = (mt_rand() % 2) == 0 ? 1 : -1;
    }

    // Initialize weights to zero
    for (int j = 0; j < num_features; ++j) {
        weights[j] = 0.0f;
    }
}

void run_computation() {
    const float learning_rate = 0.001f;
    const float lambda = 0.01f; // Regularization parameter

    // Stochastic Gradient Descent (SGD) for a linear SVM with Hinge Loss
    for (int i = 0; i < num_iterations; ++i) {
        // Select a random sample
        int sample_idx = mt_rand() % num_samples;
        float* sample_features = &data_points[sample_idx * num_features];
        int label = labels[sample_idx];

        // Calculate dot product of weights and sample features: w · x
        float dot_product = 0.0f;
        for (int j = 0; j < num_features; ++j) {
            dot_product += weights[j] * sample_features[j];
        }

        // First, apply regularization decay to all weights (gradient of regularization term)
        for (int j = 0; j < num_features; ++j) {
            weights[j] -= learning_rate * lambda * weights[j];
        }

        // Then, apply gradient update from hinge loss if misclassified or inside the margin
        // Condition for update: y * (w · x) < 1
        if (label * dot_product < 1.0f) {
            for (int j = 0; j < num_features; ++j) {
                weights[j] += learning_rate * label * sample_features[j];
            }
        }
    }

    // Calculate a final result to prevent dead code elimination
    final_result = 0.0f;
    for (int j = 0; j < num_features; ++j) {
        final_result += weights[j];
    }
}

void cleanup() {
    free(data_points);
    free(labels);
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
