#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

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
int num_samples;
int num_features;
double gamma_param;
double cost_param;
int max_iter;

// Data structures
double **X;          // Input data [num_samples][num_features]
double *y;           // Target values [num_samples]
double **K;          // Kernel matrix [num_samples][num_samples]
double *alpha;       // Lagrange multipliers [num_samples]
double *alpha_star;  // Lagrange multipliers [num_samples]
double b;            // Bias term
double final_result; // To prevent dead code elimination

// Helper to generate a random double between 0.0 and 1.0
double rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 7) {
        fprintf(stderr, "Usage: %s num_samples num_features gamma cost max_iter seed\n", argv[0]);
        exit(1);
    }

    num_samples = atoi(argv[1]);
    num_features = atoi(argv[2]);
    gamma_param = atof(argv[3]);
    cost_param = atof(argv[4]);
    max_iter = atoi(argv[5]);
    uint32_t seed = (uint32_t)atoi(argv[6]);

    mt_seed(seed);

    // Allocate memory
    X = (double **)malloc(num_samples * sizeof(double *));
    K = (double **)malloc(num_samples * sizeof(double *));
    if (!X || !K) {
         fprintf(stderr, "Memory allocation failed for matrix pointers\n"); exit(1);
    }
    for (int i = 0; i < num_samples; ++i) {
        X[i] = (double *)malloc(num_features * sizeof(double));
        K[i] = (double *)malloc(num_samples * sizeof(double));
        if (!X[i] || !K[i]) {
             fprintf(stderr, "Memory allocation failed for matrix columns\n"); exit(1);
        }
    }

    y = (double *)malloc(num_samples * sizeof(double));
    alpha = (double *)malloc(num_samples * sizeof(double));
    alpha_star = (double *)malloc(num_samples * sizeof(double));
    if (!y || !alpha || !alpha_star) {
        fprintf(stderr, "Memory allocation failed for vectors\n"); exit(1);
    }
    
    // Initialize data
    for (int i = 0; i < num_samples; ++i) {
        for (int j = 0; j < num_features; ++j) {
            X[i][j] = rand_double();
        }
        y[i] = rand_double() * 10.0; // Target values in a [-5, 5] range
        alpha[i] = 0.0;
        alpha_star[i] = 0.0;
    }
    b = 0.0;
    
    // Pre-compute the RBF kernel matrix
    for (int i = 0; i < num_samples; ++i) {
        for (int j = i; j < num_samples; ++j) {
            double dist_sq = 0.0;
            for (int k = 0; k < num_features; ++k) {
                double diff = X[i][k] - X[j][k];
                dist_sq += diff * diff;
            }
            double kernel_val = exp(-gamma_param * dist_sq);
            K[i][j] = kernel_val;
            K[j][i] = kernel_val; // Symmetric matrix
        }
    }
}

void run_computation() {
    // Simplified SVR training loop (gradient-based, not true SMO)
    const double learning_rate = 0.001;

    for (int iter = 0; iter < max_iter; ++iter) {
        for (int i = 0; i < num_samples; ++i) {
            // Calculate predicted value f(x_i)
            double f_i = 0.0;
            for (int j = 0; j < num_samples; ++j) {
                f_i += (alpha[j] - alpha_star[j]) * K[i][j];
            }
            f_i += b;

            // Calculate error
            double E_i = y[i] - f_i;

            // Update alpha and alpha_star using gradient ascent with box constraints
            alpha[i] += learning_rate * E_i;
            alpha_star[i] -= learning_rate * E_i;

            // Apply box constraints [0, C]
            alpha[i] = fmax(0.0, fmin(cost_param, alpha[i]));
            alpha_star[i] = fmax(0.0, fmin(cost_param, alpha_star[i]));

            // Update bias term b (simplified)
            b += learning_rate * E_i;
        }
    }

    // Accumulate results to prevent dead code elimination
    double result = b;
    for (int i = 0; i < num_samples; ++i) {
        result += alpha[i] - alpha_star[i];
    }
    final_result = result;
}

void cleanup() {
    for (int i = 0; i < num_samples; ++i) {
        free(X[i]);
        free(K[i]);
    }
    free(X);
    free(K);
    free(y);
    free(alpha);
    free(alpha_star);
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
