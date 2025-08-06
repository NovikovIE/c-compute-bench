#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

// Mersenne Twister (Do Not Modify)
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
// End of Mersenne Twister

// Global variables for the benchmark
int N; // Matrix size
double *A; // Input matrix, modified in-place
double *Q; // Orthogonal matrix
double *R; // Upper triangular matrix
double final_result;

// Setup function: parses arguments, allocates memory, generates data
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <matrix_size> <seed>\n", argv[0]);
        exit(1);
    }

    N = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if (N <= 0) {
        fprintf(stderr, "ERROR: matrix_size must be a positive integer.\n");
        exit(1);
    }

    // Seed the random number generator
    mt_seed(seed);

    // Allocate matrices on the heap
    // A is the input matrix, which will be modified during computation
    A = (double*)malloc(N * N * sizeof(double));
    // Q is the resulting orthogonal matrix
    Q = (double*)malloc(N * N * sizeof(double));
    // R is the resulting upper triangular matrix. calloc initializes to zero.
    R = (double*)calloc(N * N, sizeof(double));

    if (!A || !Q || !R) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }
    
    // Fill the input matrix A with random values between 0.0 and 1.0
    for (int i = 0; i < N * N; ++i) {
        A[i] = (double)mt_rand() / (double)UINT32_MAX;
    }
}

// Computation function: performs Modified Gram-Schmidt QR factorization
void run_computation() {
    // The algorithm decomposes A into Q*R, where Q is orthogonal and R is upper triangular.
    // This implementation is the Modified Gram-Schmidt process, which is more numerically stable.
    // The matrix A is consumed in the process.

    for (int k = 0; k < N; ++k) {
        // 1. Calculate the L2-norm of the k-th column vector of the working matrix A
        double norm_val = 0.0;
        for (int i = 0; i < N; ++i) {
            norm_val += A[i * N + k] * A[i * N + k];
        }
        R[k * N + k] = sqrt(norm_val);

        // 2. Compute the k-th orthogonal vector (k-th column of Q)
        // This is done by normalizing the k-th column of the current A.
        double inv_norm = 1.0 / R[k * N + k];
        for (int i = 0; i < N; ++i) {
            Q[i * N + k] = A[i * N + k] * inv_norm;
        }

        // 3. Update the remaining columns of A by subtracting the projection onto the new orthogonal vector
        for (int j = k + 1; j < N; ++j) {
            // R[k][j] is the dot product of the current j-th column of A and the k-th column of Q
            double dot_product = 0.0;
            for (int i = 0; i < N; ++i) {
                dot_product += Q[i * N + k] * A[i * N + j];
            }
            R[k * N + j] = dot_product;
            
            // Subtract the projection from the j-th column of A
            for (int i = 0; i < N; ++i) {
                A[i * N + j] -= R[k * N + j] * Q[i * N + k];
            }
        }
    }
    
    // Calculate a final result to prevent dead code elimination.
    // The trace of R (sum of diagonal elements) is a good candidate.
    final_result = 0.0;
    for (int i = 0; i < N; ++i) {
        final_result += R[i * N + i];
    }
}

// Cleanup function: frees all allocated memory
void cleanup() {
    free(A);
    free(Q);
    free(R);
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
