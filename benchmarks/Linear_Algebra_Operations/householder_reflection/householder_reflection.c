#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>
#include <limits.h>

// --- Mersenne Twister (DO NOT MODIFY) ---
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
// --- End Mersenne Twister ---

// --- Benchmark Globals ---
int N; // vector_length, matrix dimension
const int NUM_REFLECTIONS = 200; // Hardcoded number of reflections

// Input data
double *matrix_A;     // N x N matrix
double *vectors_v;    // NUM_REFLECTIONS x N matrix (each row is a vector)

// Workspace and result
double *result_matrix; // N x N matrix for computation
double *w_workspace;   // Workspace vector of size N

// Final result accumulator
double final_result;

// Helper to generate a random double in [0, 1]
double rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

// --- Benchmark Functions ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <vector_length> <seed>\n", argv[0]);
        exit(1);
    }

    N = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if (N <= 0) {
        fprintf(stderr, "vector_length must be a positive integer.\n");
        exit(1);
    }
    
    mt_seed(seed);

    size_t matrix_size = (size_t)N * N;
    size_t vectors_size = (size_t)NUM_REFLECTIONS * N;

    matrix_A = (double*)malloc(matrix_size * sizeof(double));
    vectors_v = (double*)malloc(vectors_size * sizeof(double));
    result_matrix = (double*)malloc(matrix_size * sizeof(double));
    w_workspace = (double*)malloc((size_t)N * sizeof(double));

    if (!matrix_A || !vectors_v || !result_matrix || !w_workspace) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    for (size_t i = 0; i < matrix_size; ++i) {
        matrix_A[i] = rand_double();
    }

    for (size_t i = 0; i < vectors_size; ++i) {
        vectors_v[i] = rand_double();
    }
}

void run_computation() {
    size_t matrix_bytes = (size_t)N * N * sizeof(double);
    memcpy(result_matrix, matrix_A, matrix_bytes);

    for (int r = 0; r < NUM_REFLECTIONS; ++r) {
        double* v = &vectors_v[r * N];

        // 1. Calculate norm of v squared: v_dot_v = v^T * v
        double v_dot_v = 0.0;
        for (int i = 0; i < N; ++i) {
            v_dot_v += v[i] * v[i];
        }
        
        if (v_dot_v == 0.0) continue; // Safety check for zero vector

        double beta = 2.0 / v_dot_v;

        // 2. Calculate w^T = v^T * result_matrix. Store in w_workspace.
        for (int j = 0; j < N; ++j) {
            double dot_product = 0.0;
            for (int i = 0; i < N; ++i) {
                dot_product += v[i] * result_matrix[i * N + j];
            }
            w_workspace[j] = dot_product;
        }

        // 3. Update result_matrix: A' = A - beta * v * w^T
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                result_matrix[i * N + j] -= beta * v[i] * w_workspace[j];
            }
        }
    }

    // Accumulate result to prevent dead code elimination
    double sum = 0.0;
    for (size_t i = 0; i < (size_t)N * N; ++i) {
        sum += result_matrix[i];
    }
    final_result = sum;
}

void cleanup() {
    free(matrix_A);
    free(vectors_v);
    free(result_matrix);
    free(w_workspace);
}

// --- Main ---
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
