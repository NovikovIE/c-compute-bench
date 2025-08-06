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
// --- End of Mersenne Twister ---

// Benchmark parameters and data structures
int matrix_size;
int max_iterations;
double **A; // The matrix A in the system Ax=b
double *x;  // The solution vector x, initialized to zeros
double *b;  // The vector b
double result_checksum; // Used to output a value and prevent dead-code elimination

// Generates a random double between -1.0 and 1.0
double random_double() {
    return 2.0 * ((double)mt_rand() / (double)UINT32_MAX) - 1.0;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <matrix_size> <max_iterations> <seed>\n", argv[0]);
        exit(1);
    }

    matrix_size = atoi(argv[1]);
    max_iterations = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    mt_seed(seed);

    A = (double **)malloc(matrix_size * sizeof(double *));
    if (A == NULL) { fprintf(stderr, "Memory allocation failed for A\n"); exit(1); }
    for (int i = 0; i < matrix_size; ++i) {
        A[i] = (double *)malloc(matrix_size * sizeof(double));
        if (A[i] == NULL) { fprintf(stderr, "Memory allocation failed for A[%d]\n", i); exit(1); }
    }

    x = (double *)malloc(matrix_size * sizeof(double));
    b = (double *)malloc(matrix_size * sizeof(double));
    if (x == NULL || b == NULL) { fprintf(stderr, "Memory allocation failed for x or b\n"); exit(1); }

    // Generate a diagonally dominant matrix to ensure convergence
    for (int i = 0; i < matrix_size; ++i) {
        double row_sum = 0.0;
        for (int j = 0; j < matrix_size; ++j) {
            if (i != j) {
                A[i][j] = random_double();
                row_sum += fabs(A[i][j]);
            }
        }
        // Set the diagonal element to be larger than the sum of absolute values of other elements in the row.
        A[i][i] = row_sum + 1.5; // Adding 1.5 for strong dominance
        
        // Initialize solution vector x to 0 and vector b to random values
        b[i] = random_double() * matrix_size;
        x[i] = 0.0;
    }
}

void run_computation() {
    for (int iter = 0; iter < max_iterations; ++iter) {
        for (int i = 0; i < matrix_size; ++i) {
            double sigma = 0.0;
            for (int j = 0; j < matrix_size; ++j) {
                if (i != j) {
                    sigma += A[i][j] * x[j];
                }
            }
            x[i] = (b[i] - sigma) / A[i][i];
        }
    }

    // Calculate a checksum of the result vector to prevent optimization
    result_checksum = 0.0;
    for (int i = 0; i < matrix_size; ++i) {
        result_checksum += x[i];
    }
}

void cleanup() {
    for (int i = 0; i < matrix_size; ++i) {
        free(A[i]);
    }
    free(A);
    free(x);
    free(b);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final result checksum to stdout
    printf("%f\n", result_checksum);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
