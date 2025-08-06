/*
 * Benchmark: Johnson's All-Pairs Shortest Path
 * 
 * Description: This program implements Johnson's algorithm to find the shortest paths
 * between all pairs of vertices in a weighted, directed graph. Johnson's algorithm is
 * particularly effective for sparse graphs and can handle negative edge weights,
 * provided there are no negative-weight cycles.
 *
 * The algorithm consists of three main stages:
 * 1. Reweighting using Bellman-Ford: A new vertex 's' is added to the graph with
 *    zero-weight edges to all other vertices. The Bellman-Ford algorithm is run from
 *    's' to compute a potential h(v) for each vertex v. This step also detects
 *    negative-weight cycles (though this implementation assumes none exist for simplicity).
 * 2. Edge Weight Transformation: The original edge weights w(u,v) are transformed to
 *    non-negative weights w'(u,v) = w(u,v) + h(u) - h(v).
 * 3. Dijkstra's Algorithm: Dijkstra's algorithm is run from each vertex in the original
 *    graph on the transformed weights to find the shortest paths.
 * 4. Final Path Calculation: The distances from Dijkstra's are converted back to the
 *    original shortest path distances.
 *
 * This implementation uses an adjacency list to represent the graph and a simple
 * array-based priority queue for Dijkstra's algorithm, resulting in a complexity of
 * O(V*E) for the Bellman-Ford stage and O(V^3) for the V runs of Dijkstra's, making
 * the total complexity O(V^3) for dense graphs or when using this simple Dijkstra. 
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <limits.h>
#include <stdbool.h>

// --- Mersenne Twister (MT19937) Generator - DO NOT MODIFY ---
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
// --- End of MT19937 --- 

#define INFINITY (LLONG_MAX / 2)

// --- Data Structures ---
typedef struct Edge {
    int dest;
    int weight;
    struct Edge* next;
} Edge;

// --- Global Data ---
int num_vertices;
int num_edges;

// Adjacency list for the graph. Size is num_vertices.
Edge** adj_list;

// Bellman-Ford potential values. Size is num_vertices.
long long* h_potential;

// Final all-pairs shortest path distance matrix. Size is num_vertices x num_vertices.
long long** dist_matrix;

// Final result to prevent dead code elimination.
long long final_result_checksum = 0;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_vertices> <num_edges> <seed>\n", argv[0]);
        exit(1);
    }

    num_vertices = atoi(argv[1]);
    num_edges = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    if (num_vertices <= 0 || num_edges < 0) {
        fprintf(stderr, "Error: number of vertices and edges must be positive.\n");
        exit(1);
    }

    if (num_edges > num_vertices * (num_vertices - 1)) {
        num_edges = num_vertices * (num_vertices - 1);
    }

    mt_seed(seed);

    // Allocate memory
    adj_list = (Edge**)calloc(num_vertices, sizeof(Edge*));
    h_potential = (long long*)malloc(num_vertices * sizeof(long long));
    dist_matrix = (long long**)malloc(num_vertices * sizeof(long long*));
    for (int i = 0; i < num_vertices; ++i) {
        dist_matrix[i] = (long long*)malloc(num_vertices * sizeof(long long));
    }

    if (!adj_list || !h_potential || !dist_matrix) {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        exit(1);
    }

    // Generate graph edges
    for (int i = 0; i < num_edges; ++i) {
        int u, v;
        do {
            u = mt_rand() % num_vertices;
            v = mt_rand() % num_vertices;
        } while (u == v);

        // Weight range [-10, 100]
        int weight = (mt_rand() % 111) - 10;

        Edge* new_edge = (Edge*)malloc(sizeof(Edge));
        if (!new_edge) { fprintf(stderr, "Error: Mem alloc failed for edge.\n"); exit(1); }
        new_edge->dest = v;
        new_edge->weight = weight;
        new_edge->next = adj_list[u];
        adj_list[u] = new_edge;
    }
}

void run_computation() {
    // --- Stage 1: Bellman-Ford for reweighting ---
    // Conceptually add a new source 's' with 0-weight edges to all vertices.
    for (int i = 0; i < num_vertices; ++i) {
        h_potential[i] = INFINITY;
    }
    // This is equivalent to setting h[s]=0 and relaxing all (s,v) edges once.
    for (int i = 0; i < num_vertices; ++i) h_potential[i] = 0;

    // Relax edges V-1 times
    for (int i = 0; i < num_vertices - 1; ++i) {
        for (int u = 0; u < num_vertices; ++u) {
            for (Edge* edge = adj_list[u]; edge != NULL; edge = edge->next) {
                int v = edge->dest;
                int weight = edge->weight;
                if (h_potential[u] != INFINITY && h_potential[u] + weight < h_potential[v]) {
                    h_potential[v] = h_potential[u] + weight;
                }
            }
        }
    }
    // Note: We skip the check for negative-weight cycles for simplicity in this benchmark.

    // --- Stage 2 & 3: Run Dijkstra from each vertex ---
    long long* dijkstra_dist = (long long*)malloc(num_vertices * sizeof(long long));
    bool* visited = (bool*)malloc(num_vertices * sizeof(bool));

    if (!dijkstra_dist || !visited) {
        fprintf(stderr, "Error: Malloc failed in computation.\n");
        exit(1);
    }

    for (int start_node = 0; start_node < num_vertices; ++start_node) {
        for (int i = 0; i < num_vertices; ++i) {
            dijkstra_dist[i] = INFINITY;
            visited[i] = false;
        }
        dijkstra_dist[start_node] = 0;

        for (int count = 0; count < num_vertices; ++count) {
            long long min_dist = INFINITY;
            int u = -1;

            // Find vertex with minimum distance (simple O(V) priority queue)
            for (int i = 0; i < num_vertices; ++i) {
                if (!visited[i] && dijkstra_dist[i] < min_dist) {
                    min_dist = dijkstra_dist[i];
                    u = i;
                }
            }

            if (u == -1) break; // All remaining vertices are inaccessible

            visited[u] = true;

            // Relax edges from u
            for (Edge* edge = adj_list[u]; edge != NULL; edge = edge->next) {
                int v = edge->dest;
                // Reweight edge: w'(u,v) = w(u,v) + h(u) - h(v)
                long long reweighted_w = edge->weight + h_potential[u] - h_potential[v];
                if (!visited[v] && dijkstra_dist[u] != INFINITY && dijkstra_dist[u] + reweighted_w < dijkstra_dist[v]) {
                    dijkstra_dist[v] = dijkstra_dist[u] + reweighted_w;
                }
            }
        }

        // --- Stage 4: Convert distances back ---
        for (int i = 0; i < num_vertices; ++i) {
            if (dijkstra_dist[i] == INFINITY) {
                dist_matrix[start_node][i] = INFINITY;
            } else {
                dist_matrix[start_node][i] = dijkstra_dist[i] - h_potential[start_node] + h_potential[i];
            }
        }
    }

    free(dijkstra_dist);
    free(visited);

    // --- Final checksum for verification ---
    long long checksum = 0;
    for (int i = 0; i < num_vertices; ++i) {
        for (int j = 0; j < num_vertices; ++j) {
            if (dist_matrix[i][j] != INFINITY) {
                checksum += dist_matrix[i][j];
            }
        }
    }
    final_result_checksum = checksum;
}

void cleanup() {
    for (int i = 0; i < num_vertices; ++i) {
        Edge* current = adj_list[i];
        while (current != NULL) {
            Edge* temp = current;
            current = current->next;
            free(temp);
        }
    }
    free(adj_list);

    for (int i = 0; i < num_vertices; ++i) {
        free(dist_matrix[i]);
    }
    free(dist_matrix);

    free(h_potential);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print checksum result to stdout
    printf("%lld\n", final_result_checksum);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
