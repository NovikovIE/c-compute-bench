#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <limits.h>

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
// --- End of MT19937 ---

// --- Benchmark Data Structures ---

// Represents a weighted, directed edge in the graph.
struct Edge {
    int src, dest, weight;
};

// Represents the graph using an edge list.
struct Graph {
    int v_count; // Number of vertices
    int e_count; // Number of edges
    struct Edge* edges;
};

// Global pointers to hold benchmark data
struct Graph* g_graph = NULL;
long long* g_distances = NULL;
int g_start_vertex;

// Global variable to store the final result
long long final_result = 0;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_vertices> <num_edges> <start_vertex> <seed>\n", argv[0]);
        exit(1);
    }

    int num_vertices = atoi(argv[1]);
    int num_edges = atoi(argv[2]);
    g_start_vertex = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    if (num_vertices <= 0 || num_edges < num_vertices || g_start_vertex < 0 || g_start_vertex >= num_vertices) {
        fprintf(stderr, "Invalid parameters. Note: num_edges must be >= num_vertices.\n");
        exit(1);
    }

    mt_seed(seed);

    g_graph = (struct Graph*)malloc(sizeof(struct Graph));
    if (!g_graph) { perror("Failed to allocate graph"); exit(1); }
    g_graph->v_count = num_vertices;
    g_graph->e_count = num_edges;

    g_graph->edges = (struct Edge*)malloc(num_edges * sizeof(struct Edge));
    if (!g_graph->edges) { perror("Failed to allocate edges"); free(g_graph); exit(1); }

    g_distances = (long long*)malloc(num_vertices * sizeof(long long));
    if (!g_distances) { perror("Failed to allocate distances array"); free(g_graph->edges); free(g_graph); exit(1); }

    // Generate graph data.
    // First, create a path to ensure all vertices are reachable (a cycle).
    for (int i = 0; i < num_vertices; ++i) {
        g_graph->edges[i].src = i;
        g_graph->edges[i].dest = (i + 1) % num_vertices;
        g_graph->edges[i].weight = (mt_rand() % 100) + 1; // Positive weights [1, 100]
    }

    // Then, add the remaining random edges.
    for (int i = num_vertices; i < num_edges; ++i) {
        g_graph->edges[i].src = mt_rand() % num_vertices;
        g_graph->edges[i].dest = mt_rand() % num_vertices;
        g_graph->edges[i].weight = (mt_rand() % 100) + 1; // Positive weights [1, 100]
    }
}

void run_computation() {
    int v = g_graph->v_count;
    int e = g_graph->e_count;
    struct Edge* edges = g_graph->edges;

    // Step 1: Initialize distances
    for (int i = 0; i < v; i++) {
        g_distances[i] = LLONG_MAX;
    }
    g_distances[g_start_vertex] = 0;

    // Step 2: Relax all edges |V| - 1 times.
    for (int i = 1; i <= v - 1; i++) {
        for (int j = 0; j < e; j++) {
            int u = edges[j].src;
            int v_dest = edges[j].dest;
            int weight = edges[j].weight;
            if (g_distances[u] != LLONG_MAX && g_distances[u] + weight < g_distances[v_dest]) {
                g_distances[v_dest] = g_distances[u] + weight;
            }
        }
    }

    // Since we only use positive weights, we don't need to check for negative-weight cycles.

    // Calculate a checksum to prevent dead code elimination.
    long long sum = 0;
    for (int i = 0; i < v; i++) {
        if (g_distances[i] != LLONG_MAX) {
            sum += g_distances[i];
        }
    }
    final_result = sum;
}

void cleanup() {
    if (g_graph) {
        free(g_graph->edges);
        free(g_graph);
    }
    if (g_distances) {
        free(g_distances);
    }
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
    printf("%lld\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
