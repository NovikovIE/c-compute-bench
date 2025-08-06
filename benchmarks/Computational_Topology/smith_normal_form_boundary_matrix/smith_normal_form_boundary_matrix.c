#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- Mersenne Twister (MT19937) --- (DO NOT MODIFY)
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

// --- Benchmark Globals ---
int **matrix;
int num_rows;
int num_cols;
long long final_result;

/**
 * @brief Sets up the benchmark data.
 * 
 * Parses command line arguments for matrix dimensions and random seed.
 * Allocates a 2D integer matrix to represent a boundary matrix.
 * Populates the matrix with sparse, small-magnitude random integers.
 * @param argc Argument count.
 * @param argv Argument values: num_rows, num_cols, seed.
 */
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_rows> <num_cols> <seed>\n", argv[0]);
        exit(1);
    }

    num_rows = atoi(argv[1]);
    num_cols = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (num_rows <= 0 || num_cols <= 0) {
        fprintf(stderr, "FATAL: Matrix dimensions must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    matrix = (int **)malloc(num_rows * sizeof(int *));
    if (!matrix) {
        fprintf(stderr, "FATAL: Memory allocation failed for matrix rows.\n");
        exit(1);
    }
    for (int i = 0; i < num_rows; ++i) {
        matrix[i] = (int *)malloc(num_cols * sizeof(int));
        if (!matrix[i]) {
            fprintf(stderr, "FATAL: Memory allocation failed for matrix columns.\n");
            exit(1);
        }
    }

    // Populate with sparse data, typical for boundary matrices.
    // Approx 50% sparsity and values in [-9, 9].
    for (int i = 0; i < num_rows; ++i) {
        for (int j = 0; j < num_cols; ++j) {
            if (mt_rand() % 2 == 0) {
                matrix[i][j] = 0;
            } else {
                matrix[i][j] = (int)(mt_rand() % 19) - 9;
            }
        }
    }
}

/**
 * @brief Runs the core computation: Smith Normal Form.
 * 
 * This function computes the Smith Normal Form (SNF) of the global matrix.
 * The SNF is a diagonal matrix obtained by applying integer row and column operations.
 * The algorithm iteratively finds the smallest non-zero element, moves it to the
 * pivot position, and then uses Euclidean-style reduction to zero out the rest of
 * the pivot's row and column.
 * The final result is the sum of the absolute values of the diagonal elements.
 */
void run_computation() {
    int rank_bound = num_rows < num_cols ? num_rows : num_cols;

    for (int k = 0; k < rank_bound; ++k) {
        // 1. Find pivot: element with smallest non-zero absolute value in M[k..n-1, k..m-1]
        int pivot_r = -1, pivot_c = -1;
        for (int r = k; r < num_rows; ++r) {
            for (int c = k; c < num_cols; ++c) {
                if (matrix[r][c] != 0) {
                    if (pivot_r == -1 || abs(matrix[r][c]) < abs(matrix[pivot_r][pivot_c])) {
                        pivot_r = r;
                        pivot_c = c;
                    }
                }
            }
        }

        if (pivot_r == -1) break; // Submatrix is all zeros, done.

        // 2. Move pivot to (k, k) position
        int* temp_row = matrix[k]; matrix[k] = matrix[pivot_r]; matrix[pivot_r] = temp_row;
        for (int r = 0; r < num_rows; ++r) {
            int temp_col = matrix[r][k]; matrix[r][k] = matrix[r][pivot_c]; matrix[r][pivot_c] = temp_col;
        }

        // 3. Use Euclidean algorithm via row/col operations to make the pivot the GCD
        // of its row and column, and zero out other elements.
        int reduction_needed = 1;
        while(reduction_needed) {
            reduction_needed = 0;
            // Reduce column k
            for (int r = k + 1; r < num_rows; ++r) {
                if(matrix[r][k] == 0) continue;
                int q = matrix[r][k] / matrix[k][k];
                for (int c = k; c < num_cols; ++c) { matrix[r][c] -= q * matrix[k][c]; }
                if(matrix[r][k] != 0) { // Remainder exists, swap rows to continue reduction
                    temp_row = matrix[k]; matrix[k] = matrix[r]; matrix[r] = temp_row;
                    reduction_needed = 1;
                }
            }
            // Reduce row k
            for (int c = k + 1; c < num_cols; ++c) {
                if(matrix[k][c] == 0) continue;
                int q = matrix[k][c] / matrix[k][k];
                for (int r = k; r < num_rows; ++r) { matrix[r][c] -= q * matrix[r][k]; }
                if(matrix[k][c] != 0) { // Remainder exists, swap columns
                    for(int r = 0; r < num_rows; ++r) {
                        int temp = matrix[r][k]; matrix[r][k] = matrix[r][c]; matrix[r][c] = temp;
                    }
                    reduction_needed = 1;
                }
            }
        }
    }

    final_result = 0;
    for (int i = 0; i < rank_bound; ++i) {
        final_result += (long long)abs(matrix[i][i]);
    }
}

/**
 * @brief Frees all memory allocated in `setup_benchmark`.
 */
void cleanup() {
    for (int i = 0; i < num_rows; ++i) {
        free(matrix[i]);
    }
    free(matrix);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout
    printf("%lld\n", final_result);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
