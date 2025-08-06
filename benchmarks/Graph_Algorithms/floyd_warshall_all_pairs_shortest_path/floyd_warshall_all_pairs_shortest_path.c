#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <limits.h>

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
// --- End Mersenne Twister ---

// Benchmark parameters
static int N; // Number of vertices

// Data structures
static int *dist; // Adjacency matrix for distances, stored as a 1D array

// Result
static long long final_result;

// Use a large number to represent infinity, but avoid overflow on addition.
#define INF (INT_MAX / 2)

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_vertices> <seed>\n", argv[0]);
        exit(1);
    }

    N = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if (N <= 0) {
        fprintf(stderr, "FATAL: num_vertices must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate memory for the distance matrix (flattened 2D array)
    dist = (int *)malloc((size_t)N * N * sizeof(int));
    if (dist == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed for distance matrix.\n");
        exit(1);
    }

    // Initialize the distance matrix for a dense graph
    // dist[i][j] = weight of edge from i to j
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i == j) {
                dist[i * N + j] = 0;
            } else {
                // Generate a random weight between 1 and 100
                dist[i * N + j] = (mt_rand() % 100) + 1;
            }
        }
    }
}

void run_computation() {
    // Floyd-Warshall algorithm
    for (int k = 0; k < N; ++k) {
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                // Using 1D indexing for the 2D matrix
                // Storing intermediate values can help the compiler optimize
                int d_ik = dist[i * N + k];
                int d_kj = dist[k * N + j];
                
                // If path through k is shorter, update dist[i][j]
                if (d_ik < INF && d_kj < INF) { // Check less than INF to avoid overflow
                    int new_dist = d_ik + d_kj;
                    if (new_dist < dist[i * N + j]) {
                        dist[i * N + j] = new_dist;
                    }
                }
            }
        }
    }

    // Calculate a final result to prevent dead code elimination
    final_result = 0;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            // Sum all path lengths.
            if (dist[i * N + j] < INF) {
                final_result += dist[i * N + j];
            }
        }
    }
}

void cleanup() {
    free(dist);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final accumulated result to stdout
    printf("%lld\n", final_result);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
