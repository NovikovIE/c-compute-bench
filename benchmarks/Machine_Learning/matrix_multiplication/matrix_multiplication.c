#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- Mersenne Twister (MT19937) Generator (DO NOT MODIFY) ---
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
int M, K, N;
float *A;
float *B;
float *C;
float final_result; // To prevent dead-code elimination

// --- Benchmark Functions ---

// Setup: parse arguments, allocate memory, and initialize data
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <matrix_m_dim> <matrix_k_dim> <matrix_n_dim> <seed>\n", argv[0]);
        exit(1);
    }

    M = atoi(argv[1]);
    K = atoi(argv[2]);
    N = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    
    mt_seed(seed);

    A = (float*)malloc((size_t)M * K * sizeof(float));
    B = (float*)malloc((size_t)K * N * sizeof(float));
    C = (float*)malloc((size_t)M * N * sizeof(float));

    if (!A || !B || !C) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    for (int i = 0; i < M * K; ++i) {
        A[i] = (float)mt_rand() / (float)UINT32_MAX;
    }

    for (int i = 0; i < K * N; ++i) {
        B[i] = (float)mt_rand() / (float)UINT32_MAX;
    }

    // C is the result matrix, no need to initialize with random data
}

// Computation: perform the core matrix multiplication
void run_computation() {
    // Naive matrix multiplication: C[i,j] = sum(A[i,k] * B[k,j])
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            float sum = 0.0f;
            for (int l = 0; l < K; ++l) {
                sum += A[i * K + l] * B[l * N + j];
            }
            C[i * N + j] = sum;
        }
    }

    // Calculate a checksum to prevent dead-code elimination by the compiler
    final_result = 0.0f;
    for (int i = 0; i < M * N; ++i) {
        final_result += C[i];
    }
}

// Cleanup: free all allocated memory
void cleanup() {
    free(A);
    free(B);
    free(C);
}

// --- Main and Timing --- //
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
