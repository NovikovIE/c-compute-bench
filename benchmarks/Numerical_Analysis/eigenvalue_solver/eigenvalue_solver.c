#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- Mersenne Twister (DO NOT MODIFY) ---
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
// --- End Mersenne Twister ---

// --- Benchmark Globals ---
int MATRIX_SIZE;
int NUM_ITERATIONS;
double** A;       // The matrix
double* v;        // The eigenvector guess
double* Av;       // Temporary vector for A*v product
double dominant_eigenvalue;

// --- Utility Functions ---
double random_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

// --- Benchmark Functions ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <matrix_size> <num_iterations> <seed>\n", argv[0]);
        exit(1);
    }

    MATRIX_SIZE = atoi(argv[1]);
    NUM_ITERATIONS = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    mt_seed(seed);

    // Allocate matrix and vectors
    A = (double**)malloc(MATRIX_SIZE * sizeof(double*));
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        A[i] = (double*)malloc(MATRIX_SIZE * sizeof(double));
    }
    v = (double*)malloc(MATRIX_SIZE * sizeof(double));
    Av = (double*)malloc(MATRIX_SIZE * sizeof(double));

    // Initialize a random symmetric matrix A
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        for (int j = i; j < MATRIX_SIZE; ++j) {
            double val = random_double() * 2.0 - 1.0; // Range [-1.0, 1.0]
            A[i][j] = val;
            A[j][i] = val;
        }
    }

    // Initialize a random vector v and normalize it
    double norm = 0.0;
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        v[i] = random_double();
        norm += v[i] * v[i];
    }
    norm = sqrt(norm);
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        v[i] /= norm;
    }
}

void run_computation() {
    // Power Iteration method to find the dominant eigenvalue
    for (int iter = 0; iter < NUM_ITERATIONS; ++iter) {
        // 1. Calculate the matrix-vector product Av = A * v
        for (int i = 0; i < MATRIX_SIZE; ++i) {
            Av[i] = 0.0;
            for (int j = 0; j < MATRIX_SIZE; ++j) {
                Av[i] += A[i][j] * v[j];
            }
        }

        // 2. Calculate the norm (magnitude) of Av
        // This norm is our approximation of the dominant eigenvalue
        double norm_Av = 0.0;
        for (int i = 0; i < MATRIX_SIZE; ++i) {
            norm_Av += Av[i] * Av[i];
        }
        dominant_eigenvalue = sqrt(norm_Av);

        // 3. Normalize Av to get the new v for the next iteration
        for (int i = 0; i < MATRIX_SIZE; ++i) {
            v[i] = Av[i] / dominant_eigenvalue;
        }
    }
}

void cleanup() {
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        free(A[i]);
    }
    free(A);
    free(v);
    free(Av);
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
    printf("%f\n", dominant_eigenvalue);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
