/*
 * Benchmark: discrete_morse_theory_simplification
 * Theme: Computational Topology
 * Description: This program simulates a core step in Topological Data Analysis (TDA), 
 *              specifically the simplification of a simplicial complex using a greedy
 *              algorithm inspired by Discrete Morse Theory. A simplicial complex is a
 *              high-dimensional generalization of a graph. The algorithm identifies
 *              and pairs up simplices to simplify the complex, leaving behind a set of
 *              "critical" simplices that preserve the essential topological features (homology).
 *
 *              1. Setup: A large set of abstract simplices is generated. Each simplex has a
 *                 dimension and a random value (acting as a Morse function). A sparse Hasse diagram 
 *                 (face/coface relationship graph) is built by randomly connecting simplices 
 *                 where dim(j) = dim(i) + 1.
 *
 *              2. Computation: The simplices are processed in increasing order of their values.
 *                 A greedy strategy is used: for each unpaired simplex `s`, we find the unpaired
 *                 coface `c` with the lowest value. If such a coface exists, `s` and `c` are paired
 *                 and removed from further consideration. If no such coface exists, `s` is marked
 *                 as "critical". The final result is the total count of critical simplices.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

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

// --- Benchmark Data Structures ---

typedef struct {
    int id;
    int dimension;
    uint32_t value; // Morse function value
    int paired_with; // -1 if unpaired, otherwise id of paired simplex
} Simplex;

typedef struct {
    // Parameters
    int num_simplices;
    int max_dimension;
    uint32_t seed;

    // Data structures
    Simplex* simplices;
    int** cofaces; // Adjacency list for cofaces
    int* coface_counts;
    int* sorted_indices; // Indices of simplices, sorted by value

    // Result
    int critical_simplices_count;
} BenchmarkData;

static BenchmarkData g_data;

// A constant to control the sparsity of the complex
#define CONNECTION_ATTEMPTS 30

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_simplices> <max_dimension> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_simplices = atoi(argv[1]);
    g_data.max_dimension = atoi(argv[2]);
    g_data.seed = (uint32_t)atoi(argv[3]);

    mt_seed(g_data.seed);

    // Allocate main simplex array
    g_data.simplices = (Simplex*)malloc(g_data.num_simplices * sizeof(Simplex));
    if (!g_data.simplices) { perror("malloc simplices"); exit(1); }

    // Allocate coface structures (adjacency lists)
    g_data.cofaces = (int**)malloc(g_data.num_simplices * sizeof(int*));
    g_data.coface_counts = (int*)calloc(g_data.num_simplices, sizeof(int));
    if (!g_data.cofaces || !g_data.coface_counts) { perror("malloc coface structures"); exit(1); }

    // Initialize simplices with random dimension and value
    for (int i = 0; i < g_data.num_simplices; ++i) {
        g_data.simplices[i].id = i;
        g_data.simplices[i].dimension = mt_rand() % g_data.max_dimension;
        g_data.simplices[i].value = mt_rand();
        g_data.simplices[i].paired_with = -1; // -1 means unpaired

        // Pre-allocate space for potential cofaces to simplify setup
        g_data.cofaces[i] = (int*)malloc(CONNECTION_ATTEMPTS * sizeof(int));
        if (!g_data.cofaces[i]) { perror("malloc coface list"); exit(1); }
    }

    // Build a random simplicial complex (Hasse diagram)
    // For each simplex, randomly try to connect it to potential cofaces
    for (int i = 0; i < g_data.num_simplices; ++i) {
        for (int k = 0; k < CONNECTION_ATTEMPTS; ++k) {
            int j = mt_rand() % g_data.num_simplices;
            if (i == j) continue;

            // A is a coface of B if dim(A) = dim(B)+1 and B is a face of A
            if (g_data.simplices[j].dimension == g_data.simplices[i].dimension + 1) {
                g_data.cofaces[i][g_data.coface_counts[i]++] = j;
            }
        }
    }

    // Allocate array for sorted indices (filled during computation)
    g_data.sorted_indices = (int*)malloc(g_data.num_simplices * sizeof(int));
    if (!g_data.sorted_indices) { perror("malloc sorted_indices"); exit(1); }
}

int compare_simplices(const void* a, const void* b) {
    int idx_a = *(const int*)a;
    int idx_b = *(const int*)b;
    uint32_t val_a = g_data.simplices[idx_a].value;
    uint32_t val_b = g_data.simplices[idx_b].value;
    if (val_a < val_b) return -1;
    if (val_a > val_b) return 1;
    return 0; // Should be rare with 32-bit values
}

void run_computation() {
    g_data.critical_simplices_count = 0;

    // Prepare for sorting
    for (int i = 0; i < g_data.num_simplices; ++i) {
        g_data.sorted_indices[i] = i;
    }

    // Sort simplex indices based on their value
    qsort(g_data.sorted_indices, g_data.num_simplices, sizeof(int), compare_simplices);

    // Main greedy pairing loop
    for (int i = 0; i < g_data.num_simplices; ++i) {
        int sigma_idx = g_data.sorted_indices[i];
        Simplex* sigma = &g_data.simplices[sigma_idx];

        if (sigma->paired_with != -1) {
            continue; // Already paired
        }

        // Find the unpaired coface with the minimum value
        int best_tau_idx = -1;
        uint32_t min_tau_val = UINT32_MAX;

        for (int j = 0; j < g_data.coface_counts[sigma_idx]; ++j) {
            int tau_idx = g_data.cofaces[sigma_idx][j];
            Simplex* tau = &g_data.simplices[tau_idx];

            if (tau->paired_with == -1 && tau->value < min_tau_val) {
                best_tau_idx = tau_idx;
                min_tau_val = tau->value;
            }
        }

        if (best_tau_idx != -1) {
            // Found a pair, mark both as paired
            sigma->paired_with = best_tau_idx;
            g_data.simplices[best_tau_idx].paired_with = sigma_idx;
        } else {
            // No available pair, this simplex is critical
            g_data.critical_simplices_count++;
        }
    }
}

void cleanup() {
    for (int i = 0; i < g_data.num_simplices; ++i) {
        free(g_data.cofaces[i]);
    }
    free(g_data.cofaces);
    free(g_data.coface_counts);
    free(g_data.simplices);
    free(g_data.sorted_indices);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%d\n", g_data.critical_simplices_count);
    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}
