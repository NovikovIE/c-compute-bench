#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

// --- START: Mersenne Twister (Do Not Modify) ---
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
// --- END: Mersenne Twister ---

// Helper for random doubles in [0.0, 1.0]
double rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

// Global struct for benchmark data
typedef struct {
    int num_samples;
    int num_features;
    double perplexity; // Note: Used for arg parsing, but a fixed sigma is used in this simplified impl.
    int num_iterations;
    uint32_t seed;

    double* X;          // Input data (N x F)
    double* Y;          // Output embedding (N x 2)
    double* P;          // Symmetrized high-dim affinities (N x N)
    double* D;          // Pairwise distances in high-dim (N x N)
    double* P_cond;     // Conditional probabilities P_j|i (N x N)
    double* dY;         // Gradient of cost function w.r.t. Y (N x 2)
    double* gains;      // Gains for adaptive learning rate (N x 2)
    double* Y_update;   // Momentum update term (N x 2)
    double* Q_num;      // Numerators for low-dim affinities (N x N)
    double final_result;
} BenchmarkData;

BenchmarkData g_data;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s num_samples num_features perplexity num_iterations seed\n", argv[0]);
        exit(1);
    }

    g_data.num_samples = atoi(argv[1]);
    g_data.num_features = atoi(argv[2]);
    g_data.perplexity = atof(argv[3]);
    g_data.num_iterations = atoi(argv[4]);
    g_data.seed = (uint32_t)atoi(argv[5]);

    mt_seed(g_data.seed);

    int N = g_data.num_samples;
    int F = g_data.num_features;

    // Allocate all memory required for the computation
    g_data.X = (double*)malloc(N * F * sizeof(double));
    g_data.Y = (double*)malloc(N * 2 * sizeof(double));
    g_data.P = (double*)malloc(N * N * sizeof(double));
    g_data.D = (double*)malloc(N * N * sizeof(double));
    g_data.P_cond = (double*)malloc(N * N * sizeof(double));
    g_data.dY = (double*)malloc(N * 2 * sizeof(double));
    g_data.gains = (double*)malloc(N * 2 * sizeof(double));
    g_data.Y_update = (double*)malloc(N * 2 * sizeof(double));
    g_data.Q_num = (double*)malloc(N * N * sizeof(double));

    if (!g_data.X || !g_data.Y || !g_data.P || !g_data.D || !g_data.P_cond || !g_data.dY || !g_data.gains || !g_data.Y_update || !g_data.Q_num) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(1);
    }
    
    // Initialize input data X with random values
    for (int i = 0; i < N * F; ++i) {
        g_data.X[i] = rand_double();
    }
    
    // Initialize initial embedding Y with small random values
    for (int i = 0; i < N * 2; ++i) {
        g_data.Y[i] = rand_double() * 0.0001;
    }

    // Initialize gains to 1.0, updates to 0.0
    for (int i = 0; i < N * 2; ++i) {
        g_data.gains[i] = 1.0;
        g_data.Y_update[i] = 0.0;
    }

    g_data.final_result = 0.0;
}

void run_computation() {
    int N = g_data.num_samples;
    int F = g_data.num_features;
    double* X = g_data.X;
    double* Y = g_data.Y;
    double* P = g_data.P;
    double* D = g_data.D;
    double* P_cond = g_data.P_cond;
    double* dY = g_data.dY;
    double* gains = g_data.gains;
    double* Y_update = g_data.Y_update;
    double* Q_num = g_data.Q_num;
    
    // --- Step 1: Compute pairwise affinities P --- 
    // Compute squared Euclidean distances D between high-dimensional points
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            double dist_sq = 0.0;
            for (int k = 0; k < F; ++k) {
                double diff = X[i * F + k] - X[j * F + k];
                dist_sq += diff * diff;
            }
            D[i * N + j] = dist_sq;
            D[j * N + i] = dist_sq;
        }
    }

    // Compute conditional probabilities P_j|i with fixed sigma (simplification)
    double sigma = 1.0; 
    double two_sigma_sq = 2.0 * sigma * sigma;
    for (int i = 0; i < N; ++i) {
        double row_sum = 0.0;
        for (int j = 0; j < N; ++j) {
            if (i == j) {
                P_cond[i * N + j] = 0.0;
            } else {
                P_cond[i * N + j] = exp(-D[i * N + j] / two_sigma_sq);
            }
            row_sum += P_cond[i * N + j];
        }
        if (row_sum > 1e-12) {
            for (int j = 0; j < N; ++j) {
                P_cond[i * N + j] /= row_sum;
            }
        }
    }

    // Symmetrize to get final joint probabilities P
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            P[i * N + j] = (P_cond[i * N + j] + P_cond[j * N + i]) / (2.0 * N);
        }
    }

    // --- Step 2: Gradient descent optimization --- 
    double learning_rate = 200.0;
    double min_gain = 0.01;
    
    for (int iter = 0; iter < g_data.num_iterations; ++iter) {
        double momentum = (iter < 250) ? 0.5 : 0.8;
       
        // Compute numerators of low-dimensional affinities Q
        double q_sum = 0.0;
        for (int i = 0; i < N; ++i) {
            for (int j = i + 1; j < N; ++j) {
                double diff_y0 = Y[i * 2 + 0] - Y[j * 2 + 0];
                double diff_y1 = Y[i * 2 + 1] - Y[j * 2 + 1];
                double dist_y_sq = diff_y0 * diff_y0 + diff_y1 * diff_y1;
                double num = 1.0 / (1.0 + dist_y_sq);
                Q_num[i * N + j] = Q_num[j * N + i] = num;
                q_sum += 2.0 * num;
            }
        }
        
        // Compute gradient
        for (int i = 0; i < N; ++i) {
            double grad_x = 0.0, grad_y = 0.0;
            for (int j = 0; j < N; ++j) {
                if (i == j) continue;
                double p_ij = P[i * N + j];
                double q_ij = Q_num[i * N + j] / q_sum;
                double stiffness = 4.0 * (p_ij - q_ij) * Q_num[i * N + j];
                grad_x += stiffness * (Y[i * 2 + 0] - Y[j * 2 + 0]);
                grad_y += stiffness * (Y[i * 2 + 1] - Y[j * 2 + 1]);
            }
            dY[i * 2 + 0] = grad_x;
            dY[i * 2 + 1] = grad_y;
        }

        // Update solution
        for (int i = 0; i < N * 2; ++i) {
            int sign_dY = dY[i] > 0.0;
            int sign_update = Y_update[i] > 0.0;
            gains[i] = (sign_dY == sign_update) ? (gains[i] + 0.2) : (gains[i] * 0.8);
            if (gains[i] < min_gain) gains[i] = min_gain;
            
            Y_update[i] = momentum * Y_update[i] - learning_rate * gains[i] * dY[i];
            Y[i] += Y_update[i];
        }

        // Center the solution
        double mean_y0 = 0.0, mean_y1 = 0.0;
        for (int i = 0; i < N; ++i) {
            mean_y0 += Y[i * 2 + 0];
            mean_y1 += Y[i * 2 + 1];
        }
        mean_y0 /= N;
        mean_y1 /= N;
        for (int i = 0; i < N; ++i) {
            Y[i * 2 + 0] -= mean_y0;
            Y[i * 2 + 1] -= mean_y1;
        }
    }
    
    // Accumulate the final result to prevent dead code elimination
    double sum = 0.0;
    for (int i = 0; i < N * 2; ++i) {
        sum += Y[i];
    }
    g_data.final_result = sum;
}

void cleanup() {
    free(g_data.X);
    free(g_data.Y);
    free(g_data.P);
    free(g_data.D);
    free(g_data.P_cond);
    free(g_data.dY);
    free(g_data.gains);
    free(g_data.Y_update);
    free(g_data.Q_num);
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
