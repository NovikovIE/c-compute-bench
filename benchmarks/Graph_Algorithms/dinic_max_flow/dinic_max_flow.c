#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>

// --- Mersenne Twister (MT19937) --- Do Not Modify ---
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

// --- Global Data Structures ---

#define MIN(a, b) (((a) < (b)) ? (a) : (b))

// Represents an edge in the graph.
// We implement the residual graph by directly modifying capacities.
typedef struct {
    int to;       // Destination vertex
    int capacity; // Residual capacity of the edge
    int rev;      // Index of the reverse edge in the destination's adjacency list
} Edge;

// Global structure to hold all benchmark data
static struct {
    int num_vertices;
    int source;
    int sink;

    // Adjacency list representation
    Edge** graph;
    int* edge_counts_final; // Final degrees for each vertex
    int* current_edge_indices; // To add edges during setup

    // Data for Dinic's algorithm
    int* level; // Level graph for BFS
    int* iter;  // Iterator for DFS to avoid re-exploring dead ends

    long long result_max_flow;
} g_data;

// --- Dinic's Algorithm Helpers ---

void add_edge(int from, int to, int cap) {
    int from_idx = g_data.current_edge_indices[from];
    int to_idx = g_data.current_edge_indices[to];

    // Forward edge
    g_data.graph[from][from_idx].to = to;
    g_data.graph[from][from_idx].capacity = cap;
    g_data.graph[from][from_idx].rev = to_idx;

    // Reverse edge for residual graph
    g_data.graph[to][to_idx].to = from;
    g_data.graph[to][to_idx].capacity = 0; // Directed graph, reverse capacity is initially 0
    g_data.graph[to][to_idx].rev = from_idx;

    g_data.current_edge_indices[from]++;
    g_data.current_edge_indices[to]++;
}

// Build the level graph using BFS from the source
int bfs() {
    memset(g_data.level, -1, g_data.num_vertices * sizeof(int));
    int* queue = (int*)malloc(g_data.num_vertices * sizeof(int));
    if (!queue) { perror("Failed to allocate queue for BFS"); exit(1); }
    int head = 0, tail = 0;

    g_data.level[g_data.source] = 0;
    queue[tail++] = g_data.source;

    while (head < tail) {
        int u = queue[head++];
        for (int i = 0; i < g_data.edge_counts_final[u]; ++i) {
            Edge* e = &g_data.graph[u][i];
            if (e->capacity > 0 && g_data.level[e->to] < 0) {
                g_data.level[e->to] = g_data.level[u] + 1;
                queue[tail++] = e->to;
            }
        }
    }

    free(queue);
    return g_data.level[g_data.sink] != -1;
}

// Find an augmenting path using DFS on the level graph
int dfs(int v, int t, int f) {
    if (v == t) {
        return f;
    }
    for (int* i = &g_data.iter[v]; *i < g_data.edge_counts_final[v]; (*i)++) {
        Edge* e = &g_data.graph[v][*i];
        if (e->capacity > 0 && g_data.level[v] < g_data.level[e->to]) {
            int d = dfs(e->to, t, MIN(f, e->capacity));
            if (d > 0) {
                e->capacity -= d;
                g_data.graph[e->to][e->rev].capacity += d;
                return d;
            }
        }
    }
    return 0;
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s num_vertices num_edges source_vertex sink_vertex seed\n", argv[0]);
        exit(1);
    }

    g_data.num_vertices = atoi(argv[1]);
    int num_edges = atoi(argv[2]);
    g_data.source = atoi(argv[3]);
    g_data.sink = atoi(argv[4]);
    uint32_t seed = atoi(argv[5]);

    if (g_data.num_vertices <= 0 || num_edges < 0 || g_data.source < 0 || g_data.source >= g_data.num_vertices || g_data.sink < 0 || g_data.sink >= g_data.num_vertices || g_data.source == g_data.sink) {
        fprintf(stderr, "Invalid arguments.\n");
        exit(1);
    }

    mt_seed(seed);

    // Two-pass approach: first, calculate degrees to allocate exact memory.
    g_data.edge_counts_final = (int*)calloc(g_data.num_vertices, sizeof(int));
    if (!g_data.edge_counts_final) { perror("malloc degrees"); exit(1); }

    int* temp_edges_u = malloc(num_edges * sizeof(int));
    int* temp_edges_v = malloc(num_edges * sizeof(int));
    if (!temp_edges_u || !temp_edges_v) { perror("malloc temp_edges"); exit(1); }

    for (int i = 0; i < num_edges; ++i) {
        int u = mt_rand() % g_data.num_vertices;
        int v = mt_rand() % g_data.num_vertices;
        if (u == v) {
            i--; // Avoid self-loops
            continue;
        }
        temp_edges_u[i] = u;
        temp_edges_v[i] = v;
        g_data.edge_counts_final[u]++;
        g_data.edge_counts_final[v]++; // For reverse edge
    }

    // Allocate main graph structure based on degrees
    g_data.graph = (Edge**)malloc(g_data.num_vertices * sizeof(Edge*));
    if (!g_data.graph) { perror("malloc graph"); exit(1); }
    for (int i = 0; i < g_data.num_vertices; ++i) {
        g_data.graph[i] = (Edge*)malloc(g_data.edge_counts_final[i] * sizeof(Edge));
        if (!g_data.graph[i]) { perror("malloc adjacency list"); exit(1); }
    }

    // Second pass: fill the graph with edges
    g_data.current_edge_indices = (int*)calloc(g_data.num_vertices, sizeof(int));
    if (!g_data.current_edge_indices) { perror("malloc iterators"); exit(1); }

    for (int i = 0; i < num_edges; ++i) {
        int u = temp_edges_u[i];
        int v = temp_edges_v[i];
        int cap = (mt_rand() % 100) + 1; // Capacity from 1 to 100
        add_edge(u, v, cap);
    }

    free(temp_edges_u);
    free(temp_edges_v);
    
    // Allocate algorithm-specific memory
    g_data.level = (int*)malloc(g_data.num_vertices * sizeof(int));
    g_data.iter = (int*)malloc(g_data.num_vertices * sizeof(int));
    if (!g_data.level || !g_data.iter) {
        perror("malloc level/iter");
        exit(1);
    }
}

void run_computation() {
    g_data.result_max_flow = 0;
    while (bfs()) {
        memset(g_data.iter, 0, g_data.num_vertices * sizeof(int));
        int flow;
        while ((flow = dfs(g_data.source, g_data.sink, INT_MAX)) > 0) {
            g_data.result_max_flow += flow;
        }
    }
}

void cleanup() {
    for (int i = 0; i < g_data.num_vertices; ++i) {
        free(g_data.graph[i]);
    }
    free(g_data.graph);
    free(g_data.edge_counts_final);
    free(g_data.current_edge_indices);
    free(g_data.level);
    free(g_data.iter);
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
    printf("%lld\n", g_data.result_max_flow);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
