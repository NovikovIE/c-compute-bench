#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

// --- Mersenne Twister (MT19937) --- DO NOT MODIFY ---
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
// --- End of Mersenne Twister ---

// --- Benchmark Data and Globals ---
typedef struct {
    int num_samples;
    int num_features;
    int num_epochs;
    double alpha;
    double learning_rate;
    double *X;         // Feature matrix [num_samples x num_features]
    double *y;         // Target vector [num_samples]
    double *weights;   // Model weights [num_features]
    double bias;       // Model bias
    double final_result; // To prevent dead code elimination
} BenchmarkData;

static BenchmarkData g_data;

// --- Utility Functions ---
double rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX * 2.0 - 1.0; // Range [-1.0, 1.0]
}

// --- Benchmark Functions ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 7) {
        fprintf(stderr, "Usage: %s num_samples num_features num_epochs alpha learning_rate seed\n", argv[0]);
        exit(1);
    }

    g_data.num_samples = atoi(argv[1]);
    g_data.num_features = atoi(argv[2]);
    g_data.num_epochs = atoi(argv[3]);
    g_data.alpha = atof(argv[4]);
    g_data.learning_rate = atof(argv[5]);
    uint32_t seed = (uint32_t)strtoul(argv[6], NULL, 10);

    mt_seed(seed);

    g_data.X = (double *)malloc((size_t)g_data.num_samples * g_data.num_features * sizeof(double));
    g_data.y = (double *)malloc((size_t)g_data.num_samples * sizeof(double));
    g_data.weights = (double *)malloc((size_t)g_data.num_features * sizeof(double));

    if (!g_data.X || !g_data.y || !g_data.weights) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Initialize features and targets with random data
    for (int i = 0; i < g_data.num_samples; ++i) {
        for (int j = 0; j < g_data.num_features; ++j) {
            g_data.X[i * g_data.num_features + j] = rand_double();
        }
        g_data.y[i] = rand_double();
    }

    // Initialize model parameters
    for (int j = 0; j < g_data.num_features; ++j) {
        g_data.weights[j] = 0.0;
    }
    g_data.bias = 0.0;
    g_data.final_result = 0.0;
}

void run_computation() {
    // Stochastic Gradient Descent for Ridge Regression
    for (int e = 0; e < g_data.num_epochs; ++e) {
        for (int i = 0; i < g_data.num_samples; ++i) {
            // 1. Predict y_pred
            double y_pred = g_data.bias;
            for (int j = 0; j < g_data.num_features; ++j) {
                y_pred += g_data.weights[j] * g_data.X[i * g_data.num_features + j];
            }

            // 2. Calculate error
            double error = y_pred - g_data.y[i];

            // 3. Update bias
            g_data.bias -= g_data.learning_rate * error;

            // 4. Update weights with L2 regularization term
            for (int j = 0; j < g_data.num_features; ++j) {
                double grad_weight = error * g_data.X[i * g_data.num_features + j] + 2.0 * g_data.alpha * g_data.weights[j];
                g_data.weights[j] -= g_data.learning_rate * grad_weight;
            }
        }
    }

    // Calculate a final result to prevent dead-code elimination
    double sum_weights = 0.0;
    for (int j = 0; j < g_data.num_features; ++j) {
        sum_weights += g_data.weights[j];
    }
    g_data.final_result = sum_weights + g_data.bias;
}

void cleanup() {
    free(g_data.X);
    free(g_data.y);
    free(g_data.weights);
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
    printf("%f\n", g_data.final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}