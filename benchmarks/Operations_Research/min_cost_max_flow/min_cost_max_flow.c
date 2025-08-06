#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
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

// --- BENCHMARK DATA AND ALGORITHM ---

// Edge structure for the residual graph
typedef struct {
    int to;          
    int capacity;    
    int cost;        
    int reverse_edge; 
    int next;        
} Edge;

// Global structure to hold all benchmark data
typedef struct {
    int num_nodes;
    int num_arcs;
    int source;
    int sink;
    
    Edge *edges;
    int *head;
    int edge_count;

    // Data for SPFA (Shortest Path Faster Algorithm)
    long long *dist;
    int *parent_edge;
    int *in_queue;
    int *queue;

    // Final result
    long long total_cost;
} BenchmarkData;

static BenchmarkData *data;

// Helper to add an edge and its reverse to the graph
void add_edge(int from, int to, int cap, int cost) {
    // Forward edge
    data->edges[data->edge_count].to = to;
    data->edges[data->edge_count].capacity = cap;
    data->edges[data->edge_count].cost = cost;
    data->edges[data->edge_count].reverse_edge = data->edge_count + 1;
    data->edges[data->edge_count].next = data->head[from];
    data->head[from] = data->edge_count++;

    // Reverse edge
    data->edges[data->edge_count].to = from;
    data->edges[data->edge_count].capacity = 0; // Residual capacity is initially 0
    data->edges[data->edge_count].cost = -cost; // Cost is negative for reverse edge
    data->edges[data->edge_count].reverse_edge = data->edge_count - 1;
    data->edges[data->edge_count].next = data->head[to];
    data->head[to] = data->edge_count++;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s num_nodes num_arcs seed\n", argv[0]);
        exit(1);
    }

    data = (BenchmarkData*)malloc(sizeof(BenchmarkData));
    if (!data) {
        perror("malloc failed");
        exit(1);
    }

    data->num_nodes = atoi(argv[1]);
    data->num_arcs = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    mt_seed(seed);

    data->source = 0;
    data->sink = data->num_nodes - 1;
    data->edge_count = 0;
    data->total_cost = 0;

    // Allocate memory
    data->head = (int*)malloc(data->num_nodes * sizeof(int));
    data->edges = (Edge*)malloc(data->num_arcs * 2 * sizeof(Edge));
    data->dist = (long long*)malloc(data->num_nodes * sizeof(long long));
    data->parent_edge = (int*)malloc(data->num_nodes * sizeof(int));
    data->in_queue = (int*)malloc(data->num_nodes * sizeof(int));
    data->queue = (int*)malloc(data->num_nodes * sizeof(int));
    
    if (!data->head || !data->edges || !data->dist || !data->parent_edge || !data->in_queue || !data->queue) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    for (int i = 0; i < data->num_nodes; i++) {
        data->head[i] = -1;
    }

    // Generate random graph
    for (int i = 0; i < data->num_arcs; i++) {
        int u = mt_rand() % data->num_nodes;
        int v = mt_rand() % data->num_nodes;
        if (u == v) {
            v = (v + 1) % data->num_nodes;
        }
        int capacity = (mt_rand() % 1000) + 1;
        int cost = (mt_rand() % 100) + 1;
        add_edge(u, v, capacity, cost);
    }
}

void run_computation() {
    long long total_flow = 0;
    data->total_cost = 0;

    while (1) {
        // Find shortest path in residual graph using SPFA
        for (int i = 0; i < data->num_nodes; i++) {
            data->dist[i] = LLONG_MAX;
            data->in_queue[i] = 0;
            data->parent_edge[i] = -1;
        }

        data->dist[data->source] = 0;
        int q_head = 0, q_tail = 0;
        data->queue[q_tail++] = data->source;
        data->in_queue[data->source] = 1;

        while (q_head != q_tail) {
            int u = data->queue[q_head++];
            if (q_head == data->num_nodes) q_head = 0; // Circular queue
            data->in_queue[u] = 0;

            for (int i = data->head[u]; i != -1; i = data->edges[i].next) {
                Edge *e = &data->edges[i];
                if (e->capacity > 0 && data->dist[u] != LLONG_MAX && data->dist[u] + e->cost < data->dist[e->to]) {
                    data->dist[e->to] = data->dist[u] + e->cost;
                    data->parent_edge[e->to] = i;
                    if (!data->in_queue[e->to]) {
                        data->queue[q_tail++] = e->to;
                        if (q_tail == data->num_nodes) q_tail = 0; // Circular queue
                        data->in_queue[e->to] = 1;
                    }
                }
            }
        }

        // If no augmenting path found, we are done
        if (data->dist[data->sink] == LLONG_MAX) {
            break;
        }

        // Find path flow (bottleneck capacity)
        int path_flow = INT_MAX;
        for (int cur = data->sink; cur != data->source; ) {
            int edge_idx = data->parent_edge[cur];
            Edge *e = &data->edges[edge_idx];
            if (e->capacity < path_flow) {
                path_flow = e->capacity;
            }
            cur = data->edges[e->reverse_edge].to;
        }

        // Augment flow along the path
        for (int cur = data->sink; cur != data->source; ) {
            int edge_idx = data->parent_edge[cur];
            data->edges[edge_idx].capacity -= path_flow;
            data->edges[data->edges[edge_idx].reverse_edge].capacity += path_flow;
            cur = data->edges[data->edges[edge_idx].reverse_edge].to;
        }

        total_flow += path_flow;
        data->total_cost += (long long)path_flow * data->dist[data->sink];
    }
}

void cleanup() {
    free(data->head);
    free(data->edges);
    free(data->dist);
    free(data->parent_edge);
    free(data->in_queue);
    free(data->queue);
    free(data);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%lld\n", data->total_cost);
    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}
