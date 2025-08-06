#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

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

// --- Benchmark Specific Code ---

typedef struct {
    int num_vertices;
    int source;
    int sink;
    int final_max_flow;

    // Data for Push-Relabel
    long long *excess_flow;
    int *height;
    int **capacity;
} BenchmarkData;

BenchmarkData g_data;

// Helper for min of long long and int, returns long long
long long min_val(long long a, int b) {
    return a < (long long)b ? a : (long long)b;
}

void setup_benchmark(int argc, char** argv) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <num_vertices> <num_edges> <source_vertex> <sink_vertex> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_vertices = atoi(argv[1]);
    int num_edges = atoi(argv[2]);
    g_data.source = atoi(argv[3]);
    g_data.sink = atoi(argv[4]);
    uint32_t seed = atoi(argv[5]);

    if (g_data.num_vertices <= 0 || num_edges < 0 || g_data.source < 0 || g_data.sink < 0 ||
        g_data.source >= g_data.num_vertices || g_data.sink >= g_data.num_vertices || g_data.source == g_data.sink) {
        fprintf(stderr, "Invalid arguments.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate memory for graph and algorithm data
    g_data.capacity = (int **)malloc(g_data.num_vertices * sizeof(int *));
    for (int i = 0; i < g_data.num_vertices; i++) {
        g_data.capacity[i] = (int *)calloc(g_data.num_vertices, sizeof(int));
    }
    g_data.height = (int *)calloc(g_data.num_vertices, sizeof(int));
    g_data.excess_flow = (long long *)calloc(g_data.num_vertices, sizeof(long long));
    
    // Generate graph edges with random capacities
    int max_capacity = 1000;
    for (int i = 0; i < num_edges; ++i) {
        int u, v;
        do {
            u = mt_rand() % g_data.num_vertices;
            v = mt_rand() % g_data.num_vertices;
        } while (u == v); // No self-loops

        g_data.capacity[u][v] += (mt_rand() % max_capacity) + 1;
    }

    g_data.final_max_flow = 0;
}

void cleanup() {
    for (int i = 0; i < g_data.num_vertices; i++) {
        free(g_data.capacity[i]);
    }
    free(g_data.capacity);
    free(g_data.height);
    free(g_data.excess_flow);
}

void push(int u, int v) {
    long long send = min_val(g_data.excess_flow[u], g_data.capacity[u][v]);
    g_data.capacity[u][v] -= send;
    g_data.capacity[v][u] += send;
    g_data.excess_flow[u] -= send;
    g_data.excess_flow[v] += send;
}

void relabel(int u) {
    int min_h = 2 * g_data.num_vertices; // Effectively infinity
    for (int v = 0; v < g_data.num_vertices; v++) {
        if (g_data.capacity[u][v] > 0) {
            if (g_data.height[v] < min_h) {
                min_h = g_data.height[v];
            }
        }
    }
    if (min_h < 2 * g_data.num_vertices) {
        g_data.height[u] = min_h + 1;
    }
}

void run_computation() {
    // Phase 1: Initialization
    g_data.height[g_data.source] = g_data.num_vertices;
    for (int v = 0; v < g_data.num_vertices; v++) {
        if (v == g_data.source) continue;
        if (g_data.capacity[g_data.source][v] > 0) {
            long long flow = g_data.capacity[g_data.source][v];
            push(g_data.source, v); // Using push for consistency
            g_data.excess_flow[g_data.source] += flow; // Correct excess at source
        }
    }

    // Phase 2: Main loop
    int active_node_found;
    do {
        active_node_found = 0;
        for (int u = 0; u < g_data.num_vertices; u++) {
            if (u == g_data.source || u == g_data.sink) continue;

            if (g_data.excess_flow[u] > 0) {
                active_node_found = 1;

                // Try to push flow until excess is gone or no valid neighbor exists
                while(g_data.excess_flow[u] > 0) {
                    int pushed = 0;
                    for (int v = 0; v < g_data.num_vertices; v++) {
                        if (g_data.capacity[u][v] > 0 && g_data.height[u] == g_data.height[v] + 1) {
                            push(u, v);
                            pushed = 1;
                            if (g_data.excess_flow[u] == 0) break;
                        }
                    }
                    if (!pushed) {
                        relabel(u);
                        break; // After relabeling, re-evaluate push options from the start
                    }
                }
            }
        }
    } while (active_node_found);

    g_data.final_max_flow = (int)g_data.excess_flow[g_data.sink];
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    int final_result = g_data.final_max_flow;

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%d\n", final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
