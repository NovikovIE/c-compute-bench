#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <stdbool.h>

// --- Mersenne Twister (MT19937) Generator ---
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

// --- Benchmark Data Structures ---

// This benchmark enumerates matchings in a bipartite graph.
// A matching is a set of edges without common vertices. This implementation
// finds the number of matchings that cover every vertex in the first partition (U).
// If |U| = |V|, this corresponds to enumerating perfect matchings.

typedef struct {
    int n1; // num_vertices_partition1 (U)
    int n2; // num_vertices_partition2 (V)
    double density;
    bool** adjacency_matrix; // Adjacency matrix from U to V
    unsigned long long result; // Final count of matchings
} BenchmarkData;

static BenchmarkData data;

// --- Recursive Helper for Computation ---
// This function counts the number of matchings from vertices u...n1-1 of partition 1
// to available vertices in partition 2, using backtracking.
unsigned long long count_matchings_recursive(int u, bool* matched_v) {
    // Base case: if we have successfully matched all vertices from partition 1 (U),
    // we have found one valid matching.
    if (u == data.n1) {
        return 1;
    }

    unsigned long long count = 0;
    // Iterate through all vertices 'v' in partition 2 (V)
    for (int v = 0; v < data.n2; ++v) {
        // If there is an edge from u to v and v is not yet matched
        if (data.adjacency_matrix[u][v] && !matched_v[v]) {
            // Match u with v
            matched_v[v] = true;
            
            // Recur for the next vertex in partition U
            count += count_matchings_recursive(u + 1, matched_v);
            
            // Backtrack: unmatch v to explore other possibilities
            matched_v[v] = false;
        }
    }
    return count;
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    // 1. Argument parsing
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_vertices_partition1> <num_vertices_partition2> <graph_density> <seed>\n", argv[0]);
        exit(1);
    }
    data.n1 = atoi(argv[1]);
    data.n2 = atoi(argv[2]);
    data.density = atof(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    if (data.n1 <= 0 || data.n2 <= 0 || data.density < 0.0 || data.density > 1.0 || data.n1 > 32 || data.n2 > 32) {
        fprintf(stderr, "Invalid or too large arguments.\n");
        exit(1);
    }

    // 2. Seed the random number generator
    mt_seed(seed);

    // 3. Allocate memory for the adjacency matrix
    data.adjacency_matrix = (bool**)malloc(data.n1 * sizeof(bool*));
    if (!data.adjacency_matrix) {
        perror("malloc for adjacency_matrix rows failed");
        exit(1);
    }
    for (int i = 0; i < data.n1; i++) {
        data.adjacency_matrix[i] = (bool*)malloc(data.n2 * sizeof(bool));
        if (!data.adjacency_matrix[i]) {
            perror("malloc for adjacency_matrix columns failed");
            for (int k = 0; k < i; k++) {
                free(data.adjacency_matrix[k]);
            }
            free(data.adjacency_matrix);
            exit(1);
        }
    }

    // 4. Generate graph data based on density
    uint32_t density_threshold = (uint32_t)(data.density * (double)UINT32_MAX);
    for (int i = 0; i < data.n1; i++) {
        for (int j = 0; j < data.n2; j++) {
            data.adjacency_matrix[i][j] = (mt_rand() < density_threshold);
        }
    }

    // 5. Initialize result
    data.result = 0;
}

void run_computation() {
    // A matching that covers partition 1 requires n1 <= n2.
    if (data.n1 > data.n2) {
        data.result = 0;
        return;
    }
    
    // Create a boolean array to track matched vertices in partition 2.
    // calloc initializes the memory to zero (i.e., all 'false').
    bool* matched_v = (bool*)calloc(data.n2, sizeof(bool));
    if (!matched_v) {
        perror("calloc for matched_v failed");
        exit(1);
    }

    // Start the recursive enumeration from the first vertex of partition 1 (U).
    data.result = count_matchings_recursive(0, matched_v);

    free(matched_v);
}

void cleanup() {
    if (data.adjacency_matrix) {
        for (int i = 0; i < data.n1; i++) {
            free(data.adjacency_matrix[i]);
        }
        free(data.adjacency_matrix);
        data.adjacency_matrix = NULL;
    }
}

// --- Main --- 
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%llu\n", data.result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
