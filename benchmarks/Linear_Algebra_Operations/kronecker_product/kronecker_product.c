/**
 * @file kronecker_product.c
 * @brief Benchmark for computing the Kronecker product of two dense matrices.
 *
 * This program calculates C = A ⊗ B, where A and B are square matrices.
 * The Kronecker product is a fundamental operation in linear algebra, particularly
 * in areas involving tensor products and quantum mechanics.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- BEGIN MERSENNE TWISTER (DO NOT MODIFY) ---
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
// --- END MERSENNE TWISTER ---

// Benchmark-specific global data
int M_A_SIZE, M_B_SIZE;
int M_C_DIM;

double *matrix_a;
double *matrix_b;
double *matrix_c;

double final_result;

/**
 * @brief Parses arguments, allocates memory, and initializes matrices.
 */
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <matrix_a_size> <matrix_b_size> <seed>\n", argv[0]);
        exit(1);
    }

    M_A_SIZE = atoi(argv[1]);
    M_B_SIZE = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (M_A_SIZE <= 0 || M_B_SIZE <= 0) {
        fprintf(stderr, "FATAL: Matrix sizes must be positive integers.\n");
        exit(1);
    }

    mt_seed(seed);

    M_C_DIM = M_A_SIZE * M_B_SIZE;
    size_t size_a = (size_t)M_A_SIZE * M_A_SIZE;
    size_t size_b = (size_t)M_B_SIZE * M_B_SIZE;
    size_t size_c = (size_t)M_C_DIM * M_C_DIM;

    matrix_a = (double*)malloc(size_a * sizeof(double));
    matrix_b = (double*)malloc(size_b * sizeof(double));
    matrix_c = (double*)malloc(size_c * sizeof(double));

    if (!matrix_a || !matrix_b || !matrix_c) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        free(matrix_a);
        free(matrix_b);
        free(matrix_c);
        exit(1);
    }

    // Initialize input matrices with random values between 0.0 and 1.0
    for (size_t i = 0; i < size_a; ++i) {
        matrix_a[i] = (double)mt_rand() / (double)UINT32_MAX;
    }

    for (size_t i = 0; i < size_b; ++i) {
        matrix_b[i] = (double)mt_rand() / (double)UINT32_MAX;
    }
}

/**
 * @brief Computes the Kronecker product C = A ⊗ B.
 */
void run_computation() {
    // A is M_A_SIZE x M_A_SIZE, B is M_B_SIZE x M_B_SIZE
    // C is (M_A_SIZE*M_B_SIZE) x (M_A_SIZE*M_B_SIZE)
    for (int p = 0; p < M_A_SIZE; p++) {
        for (int r = 0; r < M_A_SIZE; r++) {
            // Cache the value A[p,r] as it's used for an entire block in C
            double a_pr = matrix_a[p * M_A_SIZE + r];
            for (int q = 0; q < M_B_SIZE; q++) {
                for (int s = 0; s < M_B_SIZE; s++) {
                    // C(p*M_B_SIZE + q, r*M_B_SIZE + s) = A(p, r) * B(q, s)
                    int c_row = p * M_B_SIZE + q;
                    int c_col = r * M_B_SIZE + s;
                    matrix_c[c_row * M_C_DIM + c_col] = a_pr * matrix_b[q * M_B_SIZE + s];
                }
            }
        }
    }

    // Accumulate a result to prevent dead code elimination
    double sum = 0.0;
    size_t c_total_size = (size_t)M_C_DIM * M_C_DIM;
    for (size_t i = 0; i < c_total_size; i++) {
        sum += matrix_c[i];
    }
    final_result = sum;
}

/**
 * @brief Frees all heap-allocated memory.
 */
void cleanup() {
    free(matrix_a);
    free(matrix_b);
    free(matrix_c);
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
    printf("%f\n", final_result);

    // Print timing information to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
