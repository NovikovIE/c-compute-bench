#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- Mersenne Twister (MT19937) Generator ---
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

// --- Benchmark Globals ---
int N; // matrix_size
double **A; // Matrix A of the system Ax=b
double *b;  // Vector b of the system Ax=b
double *x;  // Solution vector x
double final_result = 0.0;

// --- Benchmark Functions ---

// Helper to generate a random double between -10.0 and 10.0
double rand_double() {
    return ((double)mt_rand() / (double)UINT32_MAX) * 20.0 - 10.0;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <matrix_size> <seed>\n", argv[0]);
        exit(1);
    }

    N = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if (N <= 0) {
        fprintf(stderr, "Error: matrix_size must be a positive integer.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate matrix A, and vectors b and x
    A = (double **)malloc(N * sizeof(double *));
    if (A == NULL) { fprintf(stderr, "Memory allocation failed for A\n"); exit(1); }
    for (int i = 0; i < N; i++) {
        A[i] = (double *)malloc(N * sizeof(double));
        if (A[i] == NULL) { fprintf(stderr, "Memory allocation failed for A[%d]\n", i); exit(1); }
    }
    b = (double *)malloc(N * sizeof(double));
    if (b == NULL) { fprintf(stderr, "Memory allocation failed for b\n"); exit(1); }
    x = (double *)malloc(N * sizeof(double));
    if (x == NULL) { fprintf(stderr, "Memory allocation failed for x\n"); exit(1); }

    // Populate matrix A and vector b with random values
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A[i][j] = rand_double();
        }
        b[i] = rand_double();
    }
}

void run_computation() {
    // Solve the linear system Ax = b using Gaussian elimination with partial pivoting

    // Forward Elimination
    for (int k = 0; k < N - 1; k++) {
        // Find pivot: row with largest element in column k (at or below k)
        int max_idx = k;
        double max_val = fabs(A[k][k]);
        for (int i = k + 1; i < N; i++) {
            if (fabs(A[i][k]) > max_val) {
                max_val = fabs(A[i][k]);
                max_idx = i;
            }
        }

        // Swap pivot row with current row (k)
        if (max_idx != k) {
            double *temp_row = A[k];
            A[k] = A[max_idx];
            A[max_idx] = temp_row;

            double temp_b = b[k];
            b[k] = b[max_idx];
            b[max_idx] = temp_b;
        }
        
        // If matrix is singular (or near-singular), pivot is zero
        if (A[k][k] == 0.0) {
            // For this benchmark, we continue, but in a real scenario, this would be an error.
            continue;
        }

        // Elimination loop
        for (int i = k + 1; i < N; i++) {
            double factor = A[i][k] / A[k][k];
            for (int j = k; j < N; j++) {
                A[i][j] -= factor * A[k][j];
            }
            b[i] -= factor * b[k];
        }
    }

    // Back Substitution
    for (int i = N - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < N; j++) {
            sum += A[i][j] * x[j];
        }
        if (A[i][i] == 0.0) { 
            x[i] = 0.0; // Avoid division by zero, result will be Na/Inf
        } else {
            x[i] = (b[i] - sum) / A[i][i];
        }
    }

    // Accumulate the results to prevent dead code elimination
    for (int i = 0; i < N; i++) {
        final_result += x[i];
    }
}

void cleanup() {
    for (int i = 0; i < N; i++) {
        free(A[i]);
    }
    free(A);
    free(b);
    free(x);
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
    printf("%f\n", final_result);

    // Print the execution time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
