#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <stdbool.h>

// Benchmark: hopcroft_karp_bipartite_matching
// This program computes the maximum cardinality matching in a bipartite graph
// using the Hopcroft-Karp algorithm.

// Defines for the Hopcroft-Karp algorithm
#define NIL 0
#define INF 0x7FFFFFFF

// --- Mersenne Twister (MT19937) a PRNG --- (DO NOT MODIFY)
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
// --- End of Mersenne Twister ---

// Global variables for the benchmark data
int num_u, num_v;    // Number of vertices in sets U (A) and V (B)
int **adj;           // Adjacency list for vertices in U
int *adj_counts;     // Number of neighbors for each vertex in U

int *pair_u, *pair_v; // Matching pairs for vertices
int *dist;           // Distance array for BFS layers
int *queue;          // Queue for BFS

int matching_size; // Final result

// Forward declarations for Hopcroft-Karp helpers
bool bfs();
bool dfs(int u);

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_vertices_set_a> <num_vertices_set_b> <num_edges> <seed>\n", argv[0]);
        exit(1);
    }

    num_u = atoi(argv[1]);
    num_v = atoi(argv[2]);
    int num_edges = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    
    mt_seed(seed);

    // Using 1-based indexing for vertices, so allocate N+1 space
    adj = (int **)malloc((num_u + 1) * sizeof(int *));
    adj_counts = (int *)calloc(num_u + 1, sizeof(int));

    // Temporary edge storage for efficient graph construction
    int(*temp_edges)[2] = malloc(num_edges * sizeof(*temp_edges));
    
    // Generate edges and count degrees
    for (int i = 0; i < num_edges; i++) {
        int u = (mt_rand() % num_u) + 1;
        int v = (mt_rand() % num_v) + 1;
        temp_edges[i][0] = u;
        temp_edges[i][1] = v;
        adj_counts[u]++;
    }

    // Allocate adjacency lists based on final degrees
    for (int i = 1; i <= num_u; i++) {
        adj[i] = (int *)malloc(adj_counts[i] * sizeof(int));
    }

    // Populate adjacency lists
    int *current_adj_indices = (int *)calloc(num_u + 1, sizeof(int));
    for (int i = 0; i < num_edges; i++) {
        int u = temp_edges[i][0];
        int v = temp_edges[i][1];
        adj[u][current_adj_indices[u]++] = v;
    }

    free(temp_edges);
    free(current_adj_indices);

    // Allocate memory for Hopcroft-Karp data structures
    pair_u = (int *)malloc((num_u + 1) * sizeof(int));
    pair_v = (int *)malloc((num_v + 1) * sizeof(int));
    dist = (int *)malloc((num_u + 1) * sizeof(int));
    queue = (int *)malloc((num_u + 1) * sizeof(int));
}

bool bfs() {
    int q_head = 0, q_tail = 0;

    for (int u = 1; u <= num_u; u++) {
        if (pair_u[u] == NIL) {
            dist[u] = 0;
            queue[q_tail++] = u;
        } else {
            dist[u] = INF;
        }
    }

    dist[NIL] = INF;

    while (q_head < q_tail) {
        int u = queue[q_head++];
        if (dist[u] < dist[NIL]) {
            for (int i = 0; i < adj_counts[u]; i++) {
                int v = adj[u][i];
                if (dist[pair_v[v]] == INF) {
                    dist[pair_v[v]] = dist[u] + 1;
                    queue[q_tail++] = pair_v[v];
                }
            }
        }
    }

    return (dist[NIL] != INF);
}

bool dfs(int u) {
    if (u != NIL) {
        for (int i = 0; i < adj_counts[u]; i++) {
            int v = adj[u][i];
            if (dist[pair_v[v]] == dist[u] + 1) {
                if (dfs(pair_v[v])) {
                    pair_v[v] = u;
                    pair_u[u] = v;
                    return true;
                }
            }
        }
        dist[u] = INF;
        return false;
    }
    return true;
}

void run_computation() {
    matching_size = 0;
    for (int i = 1; i <= num_u; i++) pair_u[i] = NIL;
    for (int i = 1; i <= num_v; i++) pair_v[i] = NIL;

    while (bfs()) {
        for (int u = 1; u <= num_u; u++) {
            if (pair_u[u] == NIL && dfs(u)) {
                matching_size++;
            }
        }
    }
}

void cleanup() {
    for (int i = 1; i <= num_u; i++) {
        free(adj[i]);
    }
    free(adj);
    free(adj_counts);
    free(pair_u);
    free(pair_v);
    free(dist);
    free(queue);
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
    printf("%d\n", matching_size);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
