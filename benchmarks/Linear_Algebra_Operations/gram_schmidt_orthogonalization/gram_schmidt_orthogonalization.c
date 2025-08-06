#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- Mersenne Twister (MT19937) Generator --- (DO NOT MODIFY)
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

// Benchmark parameters
int num_vectors;
int vector_dimension;

// Data structures
double** V; // Input vectors
double** U; // Output orthonormal vectors

double final_result; // To prevent dead code elimination

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_vectors> <vector_dimension> <seed>\n", argv[0]);
        exit(1);
    }

    num_vectors = atoi(argv[1]);
    vector_dimension = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    mt_seed(seed);

    V = (double**)malloc(num_vectors * sizeof(double*));
    U = (double**)malloc(num_vectors * sizeof(double*));
    if (!V || !U) {
        fprintf(stderr, "FATAL: Memory allocation failed for vector pointers.\n");
        exit(1);
    }

    for (int i = 0; i < num_vectors; i++) {
        V[i] = (double*)malloc(vector_dimension * sizeof(double));
        U[i] = (double*)malloc(vector_dimension * sizeof(double));
        if (!V[i] || !U[i]) {
            fprintf(stderr, "FATAL: Memory allocation failed for a vector.\n");
            exit(1);
        }
        for (int j = 0; j < vector_dimension; j++) {
            V[i][j] = ((double)mt_rand() / (double)UINT32_MAX) * 2.0 - 1.0;
        }
    }
}

void run_computation() {
    // Copy input vectors V to U to begin the in-place orthogonalization
    for (int i = 0; i < num_vectors; i++) {
        for (int j = 0; j < vector_dimension; j++) {
            U[i][j] = V[i][j];
        }
    }

    // Modified Gram-Schmidt (MGS) algorithm
    for (int k = 0; k < num_vectors; k++) {
        // 1. Normalize the k-th vector (U[k])
        double norm_sq = 0.0;
        for (int d = 0; d < vector_dimension; d++) {
            norm_sq += U[k][d] * U[k][d];
        }
        double norm = sqrt(norm_sq);

        if (norm > 1e-12) {
            for (int d = 0; d < vector_dimension; d++) {
                U[k][d] /= norm;
            }
        } else {
            // The vector is zero (or close to it), indicating linear dependence.
            // The resulting vector is set to zero.
            for (int d = 0; d < vector_dimension; d++) {
                U[k][d] = 0.0;
            }
        }

        // 2. Orthogonalize all subsequent vectors (j > k) against the new orthonormal vector U[k]
        for (int j = k + 1; j < num_vectors; j++) {
            double dot_product = 0.0;
            for (int d = 0; d < vector_dimension; d++) {
                dot_product += U[j][d] * U[k][d];
            }
            for (int d = 0; d < vector_dimension; d++) {
                U[j][d] -= dot_product * U[k][d];
            }
        }
    }

    // Calculate a checksum to prevent dead code elimination and provide a result.
    // Sum the first element of each resulting orthonormal vector.
    double checksum = 0.0;
    if (vector_dimension > 0) {
        for (int i = 0; i < num_vectors; i++) {
            checksum += U[i][0];
        }
    }
    final_result = checksum;
}

void cleanup() {
    if (V) {
        for (int i = 0; i < num_vectors; i++) {
            free(V[i]);
        }
        free(V);
    }
    if (U) {
        for (int i = 0; i < num_vectors; i++) {
            free(U[i]);
        }
        free(U);
    }
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

    // Print the timing information to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
