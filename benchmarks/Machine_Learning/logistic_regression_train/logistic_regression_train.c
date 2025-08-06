#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- Mersenne Twister (MT19937) Generator - DO NOT MODIFY ---
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
// --- End of MT19937 ---

// Benchmark parameters and data
typedef struct {
    int num_samples;
    int num_features;
    int num_epochs;
    double learning_rate;
    
    double* X;          // Feature matrix (flattened: num_samples * num_features)
    int* y;             // Labels (0 or 1)
    double* weights;
    double bias;
    
    double final_result; // To prevent dead code elimination
} BenchmarkData;

BenchmarkData g_data;

// Generates a random double in [-1.0, 1.0]
double rand_double() {
    return 2.0 * ((double)mt_rand() / (double)UINT32_MAX) - 1.0;
}

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s num_samples num_features num_epochs learning_rate seed\n", argv[0]);
        exit(1);
    }

    g_data.num_samples = atoi(argv[1]);
    g_data.num_features = atoi(argv[2]);
    g_data.num_epochs = atoi(argv[3]);
    g_data.learning_rate = atof(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    mt_seed(seed);

    // Allocate memory
    size_t x_size = (size_t)g_data.num_samples * g_data.num_features;
    g_data.X = (double*)malloc(x_size * sizeof(double));
    g_data.y = (int*)malloc((size_t)g_data.num_samples * sizeof(int));
    g_data.weights = (double*)malloc((size_t)g_data.num_features * sizeof(double));

    if (!g_data.X || !g_data.y || !g_data.weights) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Initialize data
    for (size_t i = 0; i < x_size; ++i) {
        g_data.X[i] = rand_double();
    }

    for (int i = 0; i < g_data.num_samples; ++i) {
        g_data.y[i] = mt_rand() % 2; // Binary labels: 0 or 1
    }

    for (int i = 0; i < g_data.num_features; ++i) {
        g_data.weights[i] = 0.0;
    }
    g_data.bias = 0.0;
    g_data.final_result = 0.0;
}

void run_computation() {
    for (int epoch = 0; epoch < g_data.num_epochs; ++epoch) {
        for (int i = 0; i < g_data.num_samples; ++i) {
            // Calculate linear combination z = w.T * x + b
            double z = g_data.bias;
            size_t sample_offset = (size_t)i * g_data.num_features;
            for (int j = 0; j < g_data.num_features; ++j) {
                z += g_data.weights[j] * g_data.X[sample_offset + j];
            }

            // Apply sigmoid function: prediction = 1 / (1 + exp(-z))
            double prediction = 1.0 / (1.0 + exp(-z));

            // Calculate error
            double error = prediction - g_data.y[i];

            // Update weights using gradient descent
            // dL/dw_j = (prediction - y) * x_j
            for (int j = 0; j < g_data.num_features; ++j) {
                g_data.weights[j] -= g_data.learning_rate * error * g_data.X[sample_offset + j];
            }

            // Update bias: dL/db = (prediction - y)
            g_data.bias -= g_data.learning_rate * error;
        }
    }

    // Accumulate a final result to prevent dead code elimination.
    double sum = g_data.bias;
    for (int j = 0; j < g_data.num_features; ++j) {
        sum += g_data.weights[j];
    }
    g_data.final_result = sum;
}

void cleanup() {
    free(g_data.X);
    free(g_data.y);
    free(g_data.weights);
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
    printf("%f\n", g_data.final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
