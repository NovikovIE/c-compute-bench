#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <limits.h>

// --- MERSENNE TWISTER (DO NOT MODIFY) ---
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

// --- BENCHMARK DATA AND GLOBALS ---
typedef struct {
    int num_vertices;
    int num_edges;
    int** adj;
    int* adj_counts;
    int* adj_capacity;
    int min_cycle_length;
} BenchmarkData;

static BenchmarkData g_data;

// --- HELPER FUNCTIONS ---
void add_edge(int u, int v) {
    if (g_data.adj_counts[u] >= g_data.adj_capacity[u]) {
        g_data.adj_capacity[u] = (g_data.adj_capacity[u] == 0) ? 4 : g_data.adj_capacity[u] * 2;
        g_data.adj[u] = realloc(g_data.adj[u], g_data.adj_capacity[u] * sizeof(int));
        if (!g_data.adj[u]) {
            fprintf(stderr, "Memory allocation failed for adjacency list.\n");
            exit(1);
        }
    }
    g_data.adj[u][g_data.adj_counts[u]++] = v;
}

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_vertices> <num_edges> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_vertices = atoi(argv[1]);
    g_data.num_edges = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (g_data.num_vertices <= 0 || g_data.num_edges < 0) {
        fprintf(stderr, "Number of vertices and edges must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    g_data.adj = calloc(g_data.num_vertices, sizeof(int*));
    g_data.adj_counts = calloc(g_data.num_vertices, sizeof(int));
    g_data.adj_capacity = calloc(g_data.num_vertices, sizeof(int));

    if (!g_data.adj || !g_data.adj_counts || !g_data.adj_capacity) {
        fprintf(stderr, "Memory allocation failed for graph structure.\n");
        exit(1);
    }

    for (int i = 0; i < g_data.num_edges; ++i) {
        int u, v;
        do {
            u = mt_rand() % g_data.num_vertices;
            v = mt_rand() % g_data.num_vertices;
        } while (u == v);

        add_edge(u, v);
        add_edge(v, u);
    }
}

void run_computation() {
    g_data.min_cycle_length = INT_MAX;

    int* dist = malloc(g_data.num_vertices * sizeof(int));
    int* parent = malloc(g_data.num_vertices * sizeof(int));
    int* queue = malloc(g_data.num_vertices * sizeof(int));

    if (!dist || !parent || !queue) {
        fprintf(stderr, "Failed to allocate BFS buffers.\n");
        exit(1);
    }

    for (int i = 0; i < g_data.num_vertices; i++) {
        // Run BFS from vertex i
        for (int j = 0; j < g_data.num_vertices; j++) {
            dist[j] = -1;
            parent[j] = -1;
        }

        int q_head = 0, q_tail = 0;
        queue[q_tail++] = i;
        dist[i] = 0;

        while (q_head < q_tail) {
            int u = queue[q_head++];

            for (int k = 0; k < g_data.adj_counts[u]; k++) {
                int v = g_data.adj[u][k];
                
                if (dist[v] == -1) {
                    dist[v] = dist[u] + 1;
                    parent[v] = u;
                    if (q_tail < g_data.num_vertices) {
                        queue[q_tail++] = v;
                    } else {
                        // Should not happen with sane graph sizes
                        fprintf(stderr, "BFS Queue overflow!\n");
                    }
                } else if (v != parent[u]) {
                    // Cycle detected
                    int current_cycle_len = dist[u] + dist[v] + 1;
                    if (current_cycle_len < g_data.min_cycle_length) {
                        g_data.min_cycle_length = current_cycle_len;
                    }
                }
            }
        }
    }

    free(dist);
    free(parent);
    free(queue);
}

void cleanup() {
    for (int i = 0; i < g_data.num_vertices; i++) {
        free(g_data.adj[i]);
    }
    free(g_data.adj);
    free(g_data.adj_counts);
    free(g_data.adj_capacity);
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
    printf("%d\n", g_data.min_cycle_length);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
