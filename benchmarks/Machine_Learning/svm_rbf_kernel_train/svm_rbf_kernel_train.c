#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- Mersenne Twister (MT19937) Generator (DO NOT MODIFY) ---
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

// --- Benchmark Globals and Data ---

// A struct to hold all benchmark data and parameters
static struct {
    int num_samples;
    int num_features;
    float gamma;
    float cost;

    // Data arrays
    float **X;      // Input features: [num_samples][num_features]
    int *y;         // Labels: [num_samples]
    float *alpha;   // Lagrange multipliers: [num_samples]

    // Final result to prevent dead code elimination
    float final_result_sum;
} g_data;

// Helper for checked memory allocation
void* checked_malloc(size_t size) {
    void* ptr = malloc(size);
    if (!ptr) {
        fprintf(stderr, "FATAL: malloc failed\n");
        exit(1);
    }
    return ptr;
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s num_samples num_features gamma cost seed\n", argv[0]);
        exit(1);
    }

    g_data.num_samples = atoi(argv[1]);
    g_data.num_features = atoi(argv[2]);
    g_data.gamma = atof(argv[3]);
    g_data.cost = atof(argv[4]);
    uint32_t seed = atoi(argv[5]);

    mt_seed(seed);

    // Allocate memory
    g_data.X = (float**)checked_malloc(g_data.num_samples * sizeof(float*));
    for (int i = 0; i < g_data.num_samples; ++i) {
        g_data.X[i] = (float*)checked_malloc(g_data.num_features * sizeof(float));
    }
    g_data.y = (int*)checked_malloc(g_data.num_samples * sizeof(int));
    g_data.alpha = (float*)checked_malloc(g_data.num_samples * sizeof(float));

    // Initialize data with random values
    for (int i = 0; i < g_data.num_samples; ++i) {
        for (int j = 0; j < g_data.num_features; ++j) {
            // Random float between 0.0 and 1.0
            g_data.X[i][j] = (float)mt_rand() / (float)UINT32_MAX;
        }
        // Random label of -1 or 1
        g_data.y[i] = (mt_rand() % 2) * 2 - 1;
        // Initialize alphas to 0
        g_data.alpha[i] = 0.0f;
    }
}

void cleanup() {
    for (int i = 0; i < g_data.num_samples; ++i) {
        free(g_data.X[i]);
    }
    free(g_data.X);
    free(g_data.y);
    free(g_data.alpha);
}

// RBF kernel: k(x_i, x_j) = exp(-gamma * ||x_i - x_j||^2)
static float rbf_kernel(int i, int j) {
    float dist_sq = 0.0f;
    for (int k = 0; k < g_data.num_features; ++k) {
        float diff = g_data.X[i][k] - g_data.X[j][k];
        dist_sq += diff * diff;
    }
    return expf(-g_data.gamma * dist_sq);
}

void run_computation() {
    // This is a simplified version of the SMO algorithm's main loop.
    // It doesn't aim for convergence but simulates the computational workload.
    const int num_iterations = 3; // Fixed number of passes over the dataset
    const float learning_rate = 0.01f;

    for (int iter = 0; iter < num_iterations; ++iter) {
        for (int i = 0; i < g_data.num_samples; ++i) {
            float score = 0.0f;
            // Calculate the score for sample i. This is the most expensive part.
            for (int j = 0; j < g_data.num_samples; ++j) {
                score += g_data.alpha[j] * g_data.y[j] * rbf_kernel(i, j);
            }

            // Simplified update rule based on the hinge loss gradient
            if (g_data.y[i] * score < 1.0f) {
                g_data.alpha[i] += learning_rate;
                if (g_data.alpha[i] > g_data.cost) {
                    g_data.alpha[i] = g_data.cost; // Clamp to cost
                }
            } else {
                 // Small decay if classified correctly
                g_data.alpha[i] -= learning_rate * 0.1f;
                if (g_data.alpha[i] < 0.0f) {
                    g_data.alpha[i] = 0.0f; // Clamp to 0
                }
            }
        }
    }

    // Accumulate a result to prevent dead code elimination of the alpha values
    float sum = 0.0f;
    for (int i = 0; i < g_data.num_samples; ++i) {
        sum += g_data.alpha[i];
    }
    g_data.final_result_sum = sum;
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final computational result to stdout
    printf("%f\n", g_data.final_result_sum);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}