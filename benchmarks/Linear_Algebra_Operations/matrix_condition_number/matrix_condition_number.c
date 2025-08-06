#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
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

// --- Benchmark Globals ---
int matrix_size;
double *A;
double *A_inv;
double final_result;

// Generate a random double in [0.0, 1.0]
double rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <matrix_size> <seed>\n", argv[0]);
        exit(1);
    }
    matrix_size = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);
    mt_seed(seed);

    if (matrix_size <= 0) {
        fprintf(stderr, "ERROR: matrix_size must be a positive integer.\n");
        exit(1);
    }

    A = (double*)malloc((size_t)matrix_size * matrix_size * sizeof(double));
    A_inv = (double*)malloc((size_t)matrix_size * matrix_size * sizeof(double));

    if (!A || !A_inv) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    for (int i = 0; i < matrix_size * matrix_size; ++i) {
        A[i] = rand_double() * 2.0 - 1.0; // Random values in [-1.0, 1.0]
    }
}

void run_computation() {
    int N = matrix_size;

    double *augmented = (double*)malloc((size_t)N * (size_t)(2 * N) * sizeof(double));
    if (!augmented) {
        final_result = -1.0; // Indicate error
        return;
    }

    // Initialize augmented matrix [A|I]
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            augmented[i * (2 * N) + j] = A[i * N + j];
            augmented[i * (2 * N) + N + j] = (i == j) ? 1.0 : 0.0;
        }
    }

    // Gauss-Jordan elimination
    for (int i = 0; i < N; ++i) {
        // Find pivot
        int pivot_row = i;
        for (int k = i + 1; k < N; ++k) {
            if (fabs(augmented[k * (2 * N) + i]) > fabs(augmented[pivot_row * (2 * N) + i])) {
                pivot_row = k;
            }
        }

        // Swap rows
        if (pivot_row != i) {
            for (int k = 0; k < 2 * N; ++k) {
                double temp = augmented[i * (2 * N) + k];
                augmented[i * (2 * N) + k] = augmented[pivot_row * (2 * N) + k];
                augmented[pivot_row * (2 * N) + k] = temp;
            }
        }

        // Check for singularity
        double pivot_val = augmented[i * (2 * N) + i];
        if (fabs(pivot_val) < 1e-12) {
            free(augmented);
            final_result = INFINITY;
            return;
        }

        // Normalize pivot row
        for (int j = i; j < 2 * N; ++j) {
            augmented[i * (2 * N) + j] /= pivot_val;
        }

        // Eliminate other rows
        for (int k = 0; k < N; ++k) {
            if (k != i) {
                double factor = augmented[k * (2 * N) + i];
                for (int j = i; j < 2 * N; ++j) {
                    augmented[k * (2 * N) + j] -= factor * augmented[i * (2 * N) + j];
                }
            }
        }
    }

    // Extract inverse matrix
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            A_inv[i * N + j] = augmented[i * (2 * N) + N + j];
        }
    }
    free(augmented);

    // Calculate infinity-norm of A
    double norm_A = 0.0;
    for (int i = 0; i < N; ++i) {
        double row_sum = 0.0;
        for (int j = 0; j < N; ++j) {
            row_sum += fabs(A[i * N + j]);
        }
        if (row_sum > norm_A) {
            norm_A = row_sum;
        }
    }

    // Calculate infinity-norm of A_inv
    double norm_A_inv = 0.0;
    for (int i = 0; i < N; ++i) {
        double row_sum = 0.0;
        for (int j = 0; j < N; ++j) {
            row_sum += fabs(A_inv[i * N + j]);
        }
        if (row_sum > norm_A_inv) {
            norm_A_inv = row_sum;
        }
    }

    final_result = norm_A * norm_A_inv;
}

void cleanup() {
    free(A);
    free(A_inv);
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
