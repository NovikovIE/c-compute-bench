#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- Global Benchmark Data ---
int g_matrix_rows;
int g_matrix_cols;
double *g_matrix_A;
double *g_vector_v;
double *g_vector_w; // Result vector: w = A * v
double g_final_result;

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

// Helper to generate a random double between 0.0 and 1.0
double random_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

// Function to set up data, parse arguments, and allocate memory
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <matrix_rows> <matrix_cols> <seed>\n", argv[0]);
        exit(1);
    }

    g_matrix_rows = atoi(argv[1]);
    g_matrix_cols = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (g_matrix_rows <= 0 || g_matrix_cols <= 0) {
        fprintf(stderr, "FATAL: Matrix dimensions must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate memory on the heap
    g_matrix_A = (double*)malloc((size_t)g_matrix_rows * g_matrix_cols * sizeof(double));
    g_vector_v = (double*)malloc((size_t)g_matrix_cols * sizeof(double));
    g_vector_w = (double*)malloc((size_t)g_matrix_rows * sizeof(double));

    if (!g_matrix_A || !g_vector_v || !g_vector_w) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }
    
    // Initialize matrix and vector with random values
    for (int i = 0; i < g_matrix_rows; ++i) {
        for (int j = 0; j < g_matrix_cols; ++j) {
            g_matrix_A[i * g_matrix_cols + j] = random_double();
        }
    }

    for (int i = 0; i < g_matrix_cols; ++i) {
        g_vector_v[i] = random_double();
    }
}

// Function to run the core computation
void run_computation() {
    // Perform matrix-vector multiplication: w = A * v
    for (int i = 0; i < g_matrix_rows; ++i) {
        double dot_product = 0.0;
        for (int j = 0; j < g_matrix_cols; ++j) {
            dot_product += g_matrix_A[i * g_matrix_cols + j] * g_vector_v[j];
        }
        g_vector_w[i] = dot_product;
    }

    // Accumulate the result to prevent dead code elimination
    double sum = 0.0;
    for (int i = 0; i < g_matrix_rows; ++i) {
        sum += g_vector_w[i];
    }
    g_final_result = sum;
}

// Function to free all allocated memory
void cleanup() {
    free(g_matrix_A);
    free(g_vector_v);
    free(g_vector_w);
}

// Main function: orchestrates the benchmark
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final accumulated result to stdout
    printf("%.6f\n", g_final_result);

    // Print the execution time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
