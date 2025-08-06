#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- Mersenne Twister (MT19937) Generator (Verbatim) ---
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

// Global variables for benchmark data
int N; // Matrix size
double **A; // Input matrix
double **L; // Lower triangular matrix
double **U; // Upper triangular matrix
double final_result = 0.0;

// Function Prototype for cleanup
void cleanup(void);

// Function to allocate a 2D matrix
double **allocate_matrix(int size) {
    double **matrix = (double **)malloc(size * sizeof(double *));
    if (!matrix) {
        fprintf(stderr, "FATAL: Memory allocation failed for matrix pointers.\n");
        exit(1);
    }
    for (int i = 0; i < size; i++) {
        matrix[i] = (double *)calloc(size, sizeof(double));
        if (!matrix[i]) {
            fprintf(stderr, "FATAL: Memory allocation failed for matrix row %d.\n", i);
            // Free previously allocated rows
            for(int j = 0; j < i; j++) free(matrix[j]);
            free(matrix);
            exit(1);
        }
    }
    return matrix;
}

// Function to free a 2D matrix
void free_matrix(double **matrix, int size) {
    if (matrix) {
        for (int i = 0; i < size; i++) {
            if (matrix[i]) {
                free(matrix[i]);
            }
        }
        free(matrix);
    }
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <matrix_size> <seed>\n", argv[0]);
        exit(1);
    }

    N = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if (N <= 0) {
        fprintf(stderr, "FATAL: matrix_size must be a positive integer.\n");
        exit(1);
    }

    mt_seed(seed);

    A = allocate_matrix(N);
    L = allocate_matrix(N);
    U = allocate_matrix(N);

    // Initialize matrix A with random values
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            // Scale to prevent very large/small numbers, improving stability
            A[i][j] = ((double)mt_rand() / (double)UINT32_MAX) * 10.0 + 1.0;
        }
    }
}

void run_computation() {
    // Doolittle algorithm for LU decomposition (A = LU)
    for (int k = 0; k < N; k++) {
        L[k][k] = 1.0; // Diagonal of L is 1

        // Calculate k-th row of U
        for (int j = k; j < N; j++) {
            double sum = 0.0;
            for (int p = 0; p < k; p++) {
                sum += L[k][p] * U[p][j];
            }
            U[k][j] = A[k][j] - sum;
        }

        // Check for singular matrix
        if (U[k][k] == 0.0) {
            fprintf(stderr, "FATAL: Matrix is singular. Division by zero.\n");
            // In a real scenario, pivoting would be necessary.
            // For this benchmark, we'll exit on failure.
            cleanup();
            exit(1);
        }

        // Calculate k-th column of L
        for (int i = k + 1; i < N; i++) {
            double sum = 0.0;
            for (int p = 0; p < k; p++) {
                sum += L[i][p] * U[p][k];
            }
            L[i][k] = (A[i][k] - sum) / U[k][k];
        }
    }

    // Calculate a checksum to prevent dead code elimination
    double checksum = 0.0;
    for (int i = 0; i < N; i++) {
        checksum += U[i][i]; // Sum of the diagonal of U
    }
    final_result = checksum;
}

void cleanup() {
    free_matrix(A, N);
    free_matrix(L, N);
    free_matrix(U, N);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final result to stdout
    printf("%f\n", final_result);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
