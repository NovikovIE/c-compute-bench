#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

//
// Mersenne Twister (MT19937) Generator
//
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

//
// Benchmark Globals
//

int num_vertices;
long num_edges;

struct Edge {
    int src;
    int dest;
};

struct Edge *edges;
int *parent;
int cycles_found;

//
// Disjoint Set Union (DSU) Helper Functions
//

// Finds the representative of the set containing element i (with path compression)
static int find_set(int i) {
    if (parent[i] == i)
        return i;
    return parent[i] = find_set(parent[i]);
}

// Merges the sets containing elements i and j
static void unite_sets(int i, int j) {
    int root_i = find_set(i);
    int root_j = find_set(j);
    if (root_i != root_j) {
        parent[root_i] = root_j;
    }
}

//
// Benchmark Functions
//

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_vertices> <num_edges> <seed>\n", argv[0]);
        exit(1);
    }

    num_vertices = atoi(argv[1]);
    num_edges = atol(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    mt_seed(seed);

    if (num_vertices <= 1) {
        fprintf(stderr, "Error: num_vertices must be greater than 1.\n");
        exit(1);
    }
    if (num_edges < 0) {
        fprintf(stderr, "Error: num_edges must be non-negative.\n");
        exit(1);
    }

    edges = (struct Edge*)malloc(num_edges * sizeof(struct Edge));
    parent = (int*)malloc(num_vertices * sizeof(int));

    if (!edges || !parent) {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        exit(1);
    }

    // Generate random edges. Ensure no self-loops.
    for (long i = 0; i < num_edges; ++i) {
        int u, v;
        do {
            u = mt_rand() % num_vertices;
            v = mt_rand() % num_vertices;
        } while (u == v);
        edges[i].src = u;
        edges[i].dest = v;
    }
}

void run_computation() {
    // Initialize DSU structure
    for (int i = 0; i < num_vertices; ++i) {
        parent[i] = i;
    }

    cycles_found = 0;

    // Process each edge and check for cycles
    for (long i = 0; i < num_edges; ++i) {
        int u = edges[i].src;
        int v = edges[i].dest;

        int root_u = find_set(u);
        int root_v = find_set(v);

        // If both vertices are already in the same set, adding this edge creates a cycle
        if (root_u == root_v) {
            cycles_found++;
        } else {
            // Otherwise, merge the sets
            unite_sets(u, v);
        }
    }
}

void cleanup() {
    free(edges);
    free(parent);
}

//
// Main
//
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%d\n", cycles_found);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
