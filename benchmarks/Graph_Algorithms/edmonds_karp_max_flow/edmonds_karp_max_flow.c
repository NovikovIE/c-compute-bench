#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <stdbool.h>
#include <limits.h>

// --- BEGIN MERSENNE TWISTER (DO NOT MODIFY) ---
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
// --- END MERSENNE TWISTER ---

// --- BENCHMARK DATA AND PARAMETERS ---
int V; // Number of vertices
int E; // Number of edges
int s; // Source vertex
int t; // Sink vertex

// Adjacency matrix to store residual graph capacities.
int **capacity;
int max_flow_result;

// --- BENCHMARK SETUP ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s num_vertices num_edges source_vertex sink_vertex seed\n", argv[0]);
        exit(1);
    }

    V = atoi(argv[1]);
    E = atoi(argv[2]);
    s = atoi(argv[3]);
    t = atoi(argv[4]);
    uint32_t seed = atoi(argv[5]);

    if (V <= 0 || E <= 0 || s < 0 || s >= V || t < 0 || t >= V || s == t) {
        fprintf(stderr, "Invalid parameters.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate and initialize capacity matrix (residual graph)
    capacity = (int **)calloc(V, sizeof(int *));
    if (capacity == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed for capacity matrix rows.\n");
        exit(1);
    }
    for (int i = 0; i < V; i++) {
        capacity[i] = (int *)calloc(V, sizeof(int));
        if (capacity[i] == NULL) {
            fprintf(stderr, "FATAL: Memory allocation failed for capacity matrix columns.\n");
            exit(1);
        }
    }

    // Generate random edges with random capacities
    for (int i = 0; i < E; ) {
        int u = mt_rand() % V;
        int v = mt_rand() % V;

        // Ensure no self-loops and the edge doesn't already exist.
        if (u != v && capacity[u][v] == 0) {
            // Capacity from 1 to 100
            capacity[u][v] = (mt_rand() % 100) + 1;
            i++;
        }
    }
}

// --- BENCHMARK COMPUTATION ---

// Helper function for Breadth-First Search (BFS)
// Returns true if a path from source 's' to sink 't' is found in the residual graph.
// Fills parent[] array to store the path.
static bool edmonds_karp_bfs(int *parent, int *queue, bool *visited) {
    for (int i = 0; i < V; i++) {
        visited[i] = false;
    }

    int q_front = 0, q_rear = 0;

    queue[q_rear++] = s;
    visited[s] = true;
    parent[s] = -1;

    while (q_front < q_rear) {
        int u = queue[q_front++];
        for (int v = 0; v < V; v++) {
            if (!visited[v] && capacity[u][v] > 0) {
                if (v == t) {
                    parent[v] = u;
                    return true;
                }
                queue[q_rear++] = v;
                parent[v] = u;
                visited[v] = true;
            }
        }
    }
    return false; // No path found
}

void run_computation() {
    int* parent = (int*)malloc(V * sizeof(int));
    int* queue = (int*)malloc(V * sizeof(int));
    bool* visited = (bool*)malloc(V * sizeof(bool));

    if (!parent || !queue || !visited) {
        fprintf(stderr, "FATAL: Memory allocation failed in run_computation.\n");
        // Cleanup is not called here, but the program will exit.
        // In a real application, memory should be handled more gracefully.
        exit(1);
    }

    max_flow_result = 0;
    
    // While there is an augmenting path from source to sink
    while (edmonds_karp_bfs(parent, queue, visited)) {
        int path_flow = INT_MAX;

        // Find the maximum flow through the path found by BFS
        for (int v = t; v != s; v = parent[v]) {
            int u = parent[v];
            if (capacity[u][v] < path_flow) {
                path_flow = capacity[u][v];
            }
        }

        // Update residual capacities of the edges and reverse edges along the path
        for (int v = t; v != s; v = parent[v]) {
            int u = parent[v];
            capacity[u][v] -= path_flow;
            capacity[v][u] += path_flow;
        }

        // Add path flow to overall flow
        max_flow_result += path_flow;
    }

    free(parent);
    free(queue);
    free(visited);
}

// --- BENCHMARK CLEANUP ---
void cleanup() {
    if (capacity) {
        for (int i = 0; i < V; i++) {
            free(capacity[i]);
        }
        free(capacity);
    }
}

// --- MAIN FUNCTION ---
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final result to stdout
    printf("%d\n", max_flow_result);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
