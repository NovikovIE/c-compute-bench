#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

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


// --- Benchmark Data Structures ---
typedef struct {
    int src, dest, weight;
} Edge;

// Global parameters
int NUM_VERTICES;
int NUM_EDGES;

// Global data structures
Edge *graph_edges;
int *parent; // For Disjoint Set Union (DSU)

// Global result
long long total_mst_weight;

// --- DSU Helper Functions ---
// Find set of an element i (with path compression)
int find_set(int i) {
    if (parent[i] == i)
        return i;
    // Path compression
    return parent[i] = find_set(parent[i]);
}

// Union of two sets x and y (simple union)
void union_sets(int x, int y) {
    int root_x = find_set(x);
    int root_y = find_set(y);
    if (root_x != root_y) {
        parent[root_x] = root_y;
    }
}

// Comparison function for qsort to sort edges by weight
int compare_edges(const void* a, const void* b) {
    Edge* edge_a = (Edge*)a;
    Edge* edge_b = (Edge*)b;
    if (edge_a->weight < edge_b->weight) return -1;
    if (edge_a->weight > edge_b->weight) return 1;
    return 0;
}


// --- Benchmark Functions ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_vertices> <num_edges> <seed>\n", argv[0]);
        exit(1);
    }
    NUM_VERTICES = atoi(argv[1]);
    NUM_EDGES = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if(NUM_VERTICES <= 0 || NUM_EDGES <= 0) {
        fprintf(stderr, "FATAL: num_vertices and num_edges must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate memory for the graph (edge list)
    graph_edges = (Edge *)malloc(NUM_EDGES * sizeof(Edge));
    if (!graph_edges) {
        fprintf(stderr, "FATAL: Memory allocation failed for graph_edges.\n");
        exit(1);
    }

    // Allocate memory for the DSU parent array
    parent = (int *)malloc(NUM_VERTICES * sizeof(int));
    if (!parent) {
        fprintf(stderr, "FATAL: Memory allocation failed for parent array.\n");
        free(graph_edges);
        exit(1);
    }

    // Generate random edges for the graph
    for (int i = 0; i < NUM_EDGES; ++i) {
        int u = mt_rand() % NUM_VERTICES;
        int v = mt_rand() % NUM_VERTICES;
        // Ensure u and v are different to avoid self-loops.
        while (u == v) {
            v = mt_rand() % NUM_VERTICES;
        }
        graph_edges[i].src = u;
        graph_edges[i].dest = v;
        // Generate a positive weight
        graph_edges[i].weight = (int)(mt_rand() % 10000) + 1;
    }
}

void run_computation() {
    // 1. Sort all the edges in non-decreasing order of their weight.
    qsort(graph_edges, NUM_EDGES, sizeof(Edge), compare_edges);

    // 2. Initialize DSU: make each vertex its own disjoint set.
    for (int i = 0; i < NUM_VERTICES; ++i) {
        parent[i] = i;
    }

    // 3. Iterate through sorted edges and build the MST
    total_mst_weight = 0;
    int mst_edge_count = 0;
    int edge_idx = 0;

    // The MST will have V-1 edges for a connected graph
    while (mst_edge_count < NUM_VERTICES - 1 && edge_idx < NUM_EDGES) {
        Edge next_edge = graph_edges[edge_idx++];

        int set_u = find_set(next_edge.src);
        int set_v = find_set(next_edge.dest);

        // If including this edge doesn't cause a cycle, include it
        if (set_u != set_v) {
            total_mst_weight += next_edge.weight;
            union_sets(set_u, set_v);
            mst_edge_count++;
        }
    }
}

void cleanup() {
    free(graph_edges);
    free(parent);
}

// --- Main Function ---
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%lld\n", total_mst_weight);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
