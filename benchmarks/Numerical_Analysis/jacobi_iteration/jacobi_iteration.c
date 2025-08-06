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

// Global variables for benchmark data
int matrix_size;
int max_iterations;
double **A;      // Coefficient matrix
double *b;       // Constant vector
double *x;       // Solution vector (current)
double *x_new;    // Solution vector (next)
double final_result; // To prevent dead code elimination

// Generates a random double between 0.0 and 1.0
double rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <matrix_size> <max_iterations> <seed>\n", argv[0]);
        exit(1);
    }

    matrix_size = atoi(argv[1]);
    max_iterations = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (matrix_size <= 0 || max_iterations <= 0) {
        fprintf(stderr, "Matrix size and max iterations must be positive integers.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate memory for matrices and vectors
    A = (double **)malloc(matrix_size * sizeof(double *));
    for (int i = 0; i < matrix_size; i++) {
        A[i] = (double *)malloc(matrix_size * sizeof(double));
    }
    b = (double *)malloc(matrix_size * sizeof(double));
    x = (double *)malloc(matrix_size * sizeof(double));
    x_new = (double *)malloc(matrix_size * sizeof(double));

    // Generate a diagonally dominant matrix A and vector b to ensure convergence
    for (int i = 0; i < matrix_size; i++) {
        double row_sum = 0.0;
        for (int j = 0; j < matrix_size; j++) {
            if (i != j) {
                A[i][j] = rand_double();
                row_sum += fabs(A[i][j]);
            }
        }
        A[i][i] = row_sum + rand_double() + 1.0; // Ensure diagonal dominance
        b[i] = rand_double() * matrix_size;
        x[i] = 0.0; // Initial guess for x is all zeros
    }
}

void run_computation() {
    for (int iter = 0; iter < max_iterations; iter++) {
        for (int i = 0; i < matrix_size; i++) {
            double sigma = 0.0;
            for (int j = 0; j < matrix_size; j++) {
                if (i != j) {
                    sigma += A[i][j] * x[j];
                }
            }
            x_new[i] = (b[i] - sigma) / A[i][i];
        }

        // Update x with the new values
        for (int i = 0; i < matrix_size; i++) {
            x[i] = x_new[i];
        }
    }

    // Calculate a final sum to report as the result and prevent optimization
    double sum = 0.0;
    for (int i = 0; i < matrix_size; i++) {
        sum += x[i];
    }
    final_result = sum;
}

void cleanup() {
    for (int i = 0; i < matrix_size; i++) {
        free(A[i]);
    }
    free(A);
    free(b);
    free(x);
    free(x_new);
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
    printf("%f\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
