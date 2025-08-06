#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>
#include <stdbool.h>

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

// --- Global Benchmark Data ---
int num_vertices;
int num_edges;

// Adjacency list representation for an undirected graph
int** adj_list; 
int* adj_list_sizes;

double final_result; // Accumulated closeness centrality values

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_vertices> <num_edges> <seed>\n", argv[0]);
        exit(1);
    }

    num_vertices = atoi(argv[1]);
    num_edges = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);
    mt_seed(seed);

    if (num_vertices <= 0 || num_edges < 0) {
        fprintf(stderr, "Number of vertices must be positive, and edges non-negative.\n");
        exit(1);
    }

    // To represent edges and degrees during construction
    int* degrees = (int*)calloc(num_vertices, sizeof(int));
    bool** edge_exists = (bool**)malloc(num_vertices * sizeof(bool*));
    for(int i = 0; i < num_vertices; ++i) {
        edge_exists[i] = (bool*)calloc(num_vertices, sizeof(bool));
    }

    int edges_added = 0;

    // 1. Create a random spanning tree to ensure the graph is connected.
    if (num_vertices > 1) {
        int* p = (int*)malloc(num_vertices * sizeof(int));
        for (int i = 0; i < num_vertices; i++) p[i] = i;
        for (int i = num_vertices - 1; i > 0; i--) { // Fisher-Yates shuffle
            int j = mt_rand() % (i + 1);
            int temp = p[i]; p[i] = p[j]; p[j] = temp;
        }
        for (int i = 0; i < num_vertices - 1; i++) {
            int u = p[i];
            int v = p[i + 1];
            edge_exists[u][v] = edge_exists[v][u] = true;
            degrees[u]++;
            degrees[v]++;
            edges_added++;
        }
        free(p);
    }

    // 2. Add remaining edges randomly, avoiding self-loops and duplicates.
    long long max_edges = (long long)num_vertices * (num_vertices - 1) / 2;
    if (num_edges > max_edges) num_edges = max_edges;

    while (edges_added < num_edges) {
        int u = mt_rand() % num_vertices;
        int v = mt_rand() % num_vertices;
        if (u == v || edge_exists[u][v]) {
            continue;
        }
        edge_exists[u][v] = edge_exists[v][u] = true;
        degrees[u]++;
        degrees[v]++;
        edges_added++;
    }

    // 3. Allocate and populate the final adjacency list.
    adj_list = (int**)malloc(num_vertices * sizeof(int*));
    adj_list_sizes = (int*)malloc(num_vertices * sizeof(int));
    int* current_adj_idx = (int*)calloc(num_vertices, sizeof(int));

    memcpy(adj_list_sizes, degrees, num_vertices * sizeof(int));
    
    for (int i = 0; i < num_vertices; i++) {
        adj_list[i] = (int*)malloc(adj_list_sizes[i] * sizeof(int));
    }

    for (int u = 0; u < num_vertices; u++) {
        for (int v = u + 1; v < num_vertices; v++) {
            if (edge_exists[u][v]) {
                adj_list[u][current_adj_idx[u]++] = v;
                adj_list[v][current_adj_idx[v]++] = u;
            }
        }
    }

    // Free temporary structures
    free(degrees);
    free(current_adj_idx);
    for (int i = 0; i < num_vertices; i++) {
        free(edge_exists[i]);
    }
    free(edge_exists);
}

void run_computation() {
    final_result = 0.0;
    if (num_vertices <= 1) return;

    // Allocate memory for BFS (Breadth-First Search)
    int* dist = (int*)malloc(num_vertices * sizeof(int));
    int* queue = (int*)malloc(num_vertices * sizeof(int));

    // Run BFS from each vertex to find all-pairs shortest paths
    for (int start_node = 0; start_node < num_vertices; start_node++) {
        // Initialize distances
        for (int i = 0; i < num_vertices; i++) dist[i] = -1;

        // BFS setup
        int head = 0, tail = 0;
        dist[start_node] = 0;
        queue[tail++] = start_node;

        while (head < tail) {
            int u = queue[head++];
            for (int i = 0; i < adj_list_sizes[u]; i++) {
                int v = adj_list[u][i];
                if (dist[v] == -1) {
                    dist[v] = dist[u] + 1;
                    queue[tail++] = v;
                }
            }
        }

        // Calculate sum of distances for the current start_node
        long long sum_of_distances = 0;
        for (int i = 0; i < num_vertices; i++) {
            if (dist[i] != -1) { // Only consider reachable nodes
                sum_of_distances += dist[i];
            }
        }

        // Calculate closeness centrality and accumulate
        if (sum_of_distances > 0) {
            // Using harmonic centrality as the computed value
            double centrality = 1.0 / (double)sum_of_distances;
            final_result += centrality;
        }
    }

    free(dist);
    free(queue);
}

void cleanup() {
    for (int i = 0; i < num_vertices; i++) {
        free(adj_list[i]);
    }
    free(adj_list);
    free(adj_list_sizes);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%f\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
