#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <limits.h>

// --- START MERSENNE TWISTER (MT19937) --- Do Not Modify ---
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
// --- END MERSENNE TWISTER ---

// --- BENCHMARK DATA AND GLOBALS ---
int N; // Number of matrices
int *p; // Array of dimensions p[0], p[1], ..., p[N]. Size N+1
long long **m; // DP table to store results of subproblems. Size N x N
long long final_result; // The minimum number of scalar multiplications

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_matrices> <seed>\n", argv[0]);
        exit(1);
    }

    N = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if (N <= 0) {
        fprintf(stderr, "Error: num_matrices must be a positive integer.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate memory for the dimensions array p
    // There are N matrices, so N+1 dimensions
    p = (int *)malloc((N + 1) * sizeof(int));
    if (p == NULL) {
        fprintf(stderr, "Failed to allocate memory for dimensions array.\n");
        exit(1);
    }

    // Generate random dimensions for matrices. Range [10, 110]
    for (int i = 0; i <= N; i++) {
        p[i] = (mt_rand() % 101) + 10;
    }

    // Allocate memory for the DP table m
    m = (long long **)malloc(N * sizeof(long long *));
    if (m == NULL) {
        fprintf(stderr, "Failed to allocate memory for DP table rows.\n");
        free(p);
        exit(1);
    }
    for (int i = 0; i < N; i++) {
        m[i] = (long long *)malloc(N * sizeof(long long));
        if (m[i] == NULL) {
            fprintf(stderr, "Failed to allocate memory for DP table columns.\n");
            // Cleanup previously allocated rows
            for (int k = 0; k < i; k++) {
                free(m[k]);
            }
            free(m);
            free(p);
            exit(1);
        }
    }
}

void run_computation() {
    // Matrix Ai has dimension p[i] x p[i+1] for i = 0..N-1
    // The DP table m[i][j] stores the minimum number of multiplications
    // needed to compute the product of matrices from A_i to A_j.

    // A single matrix requires 0 multiplications.
    for (int i = 0; i < N; i++) {
        m[i][i] = 0;
    }

    // L is the length of the matrix chain.
    for (int L = 2; L <= N; L++) {
        // i is the starting index of the chain.
        for (int i = 0; i <= N - L; i++) {
            // j is the ending index of the chain.
            int j = i + L - 1;
            m[i][j] = LLONG_MAX;
            // k is the split point.
            for (int k = i; k < j; k++) {
                // Cost = cost to compute M(i..k) + cost to compute M(k+1..j) + cost to multiply the results.
                // M(i..k) has dimensions p[i] x p[k+1]
                // M(k+1..j) has dimensions p[k+1] x p[j+1]
                long long cost = m[i][k] + m[k+1][j] + (long long)p[i] * p[k+1] * p[j+1];
                if (cost < m[i][j]) {
                    m[i][j] = cost;
                }
            }
        }
    }

    final_result = m[0][N - 1];
}

void cleanup() {
    if (m != NULL) {
        for (int i = 0; i < N; i++) {
            if (m[i] != NULL) {
                free(m[i]);
            }
        }
        free(m);
    }
    if (p != NULL) {
        free(p);
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
    printf("%lld\n", final_result);

    // Print the execution time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}