#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- BEGIN: Mersenne Twister (DO NOT MODIFY) ---
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

// --- Global Benchmark Data ---
typedef struct {
    int matrix_rows;
    int matrix_cols;
    int target_rank;
    double* matrix;
    int computed_rank;
} BenchmarkData;

BenchmarkData g_data;
const double EPSILON = 1e-9; // Tolerance for floating point comparisons

// Helper to generate a random double between 0.0 and 1.0
double rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <matrix_rows> <matrix_cols> <target_rank> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.matrix_rows = atoi(argv[1]);
    g_data.matrix_cols = atoi(argv[2]);
    g_data.target_rank = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    mt_seed(seed);

    if (g_data.target_rank > g_data.matrix_rows || g_data.target_rank > g_data.matrix_cols) {
        fprintf(stderr, "FATAL: target_rank cannot be greater than matrix dimensions.\n");
        exit(1);
    }

    // To generate a matrix with a specific rank, we create two smaller matrices, U and V,
    // and compute their product M = U * V. The rank of M will be at most the inner dimension,
    // which is the target_rank. With random values, the rank is very likely to be exactly target_rank.
    double* U = (double*)malloc((size_t)g_data.matrix_rows * g_data.target_rank * sizeof(double));
    double* V = (double*)malloc((size_t)g_data.target_rank * g_data.matrix_cols * sizeof(double));
    g_data.matrix = (double*)malloc((size_t)g_data.matrix_rows * g_data.matrix_cols * sizeof(double));

    if (!U || !V || !g_data.matrix) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        free(U); free(V); free(g_data.matrix);
        exit(1);
    }

    // Fill U and V with random values
    for (int i = 0; i < g_data.matrix_rows * g_data.target_rank; ++i) {
        U[i] = rand_double();
    }
    for (int i = 0; i < g_data.target_rank * g_data.matrix_cols; ++i) {
        V[i] = rand_double();
    }

    // Compute M = U * V
    for (int i = 0; i < g_data.matrix_rows; ++i) {
        for (int j = 0; j < g_data.matrix_cols; ++j) {
            double sum = 0.0;
            for (int k = 0; k < g_data.target_rank; ++k) {
                sum += U[i * g_data.target_rank + k] * V[k * g_data.matrix_cols + j];
            }
            g_data.matrix[i * g_data.matrix_cols + j] = sum;
        }
    }

    free(U);
    free(V);
}

void run_computation() {
    int rank = 0;
    int pivot_row = 0;
    
    // Use Gaussian elimination to convert the matrix to row-echelon form.
    // The number of non-zero rows (or equivalently, the number of pivots) in the echelon form is the rank.
    for (int col = 0; col < g_data.matrix_cols && pivot_row < g_data.matrix_rows; ++col) {
        // Find a row with a non-zero entry in the current column (the pivot)
        int p = pivot_row;
        while (p < g_data.matrix_rows && fabs(g_data.matrix[p * g_data.matrix_cols + col]) < EPSILON) {
            p++;
        }

        if (p < g_data.matrix_rows) { // A pivot was found
            // Swap the pivot row with the current pivot_row position
            if (p != pivot_row) {
                for (int k = 0; k < g_data.matrix_cols; k++) {
                    double temp = g_data.matrix[pivot_row * g_data.matrix_cols + k];
                    g_data.matrix[pivot_row * g_data.matrix_cols + k] = g_data.matrix[p * g_data.matrix_cols + k];
                    g_data.matrix[p * g_data.matrix_cols + k] = temp;
                }
            }

            // Eliminate all entries below the pivot in the current column
            for (int i = pivot_row + 1; i < g_data.matrix_rows; ++i) {
                double factor = g_data.matrix[i * g_data.matrix_cols + col] / g_data.matrix[pivot_row * g_data.matrix_cols + col];
                for (int k = col; k < g_data.matrix_cols; ++k) {
                    g_data.matrix[i * g_data.matrix_cols + k] -= factor * g_data.matrix[pivot_row * g_data.matrix_cols + k];
                }
            }
            pivot_row++;
        }
    }

    g_data.computed_rank = pivot_row;
}

void cleanup() {
    if (g_data.matrix) {
        free(g_data.matrix);
        g_data.matrix = NULL;
    }
}

// --- Main Execution Logic ---

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout
    printf("%d\n", g_data.computed_rank);

    // Print timing information to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
