#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

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

// Global variables for the benchmark
static int N_VERTICES;
static int M_EDGES;
static int** adj;
static int* out_degree;
static int* in_degree;
static int* topological_sort_result;
static long long final_result;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_vertices> <num_edges> <seed>\n", argv[0]);
        exit(1);
    }

    N_VERTICES = atoi(argv[1]);
    M_EDGES = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);
    mt_seed(seed);

    if (N_VERTICES <= 1) {
        fprintf(stderr, "Number of vertices must be greater than 1.\n");
        exit(1);
    }

    // Allocate memory for graph data
    adj = (int**)malloc(N_VERTICES * sizeof(int*));
    out_degree = (int*)calloc(N_VERTICES, sizeof(int));
    in_degree = (int*)calloc(N_VERTICES, sizeof(int));
    topological_sort_result = (int*)malloc(N_VERTICES * sizeof(int));
    
    // Step 1: Generate edges and count degrees.
    // To guarantee a DAG, we ensure that for every edge (u, v), u < v.
    int (*edge_list)[2] = malloc(M_EDGES * sizeof(*edge_list));
    if (!edge_list || !adj || !out_degree || !in_degree || !topological_sort_result) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(1);
    }
    
    for (int i = 0; i < M_EDGES; i++) {
        int u, v;
        do {
            u = mt_rand() % N_VERTICES;
            v = mt_rand() % N_VERTICES;
        } while (u == v);

        if (u > v) {
            int temp = u;
            u = v;
            v = temp;
        }
        
        edge_list[i][0] = u;
        edge_list[i][1] = v;
        out_degree[u]++;
        in_degree[v]++;
    }

    // Step 2: Allocate memory for adjacency lists based on out-degrees.
    for (int i = 0; i < N_VERTICES; i++) {
        if (out_degree[i] > 0) {
            adj[i] = (int*)malloc(out_degree[i] * sizeof(int));
        } else {
            adj[i] = NULL;
        }
    }

    // Step 3: Populate adjacency lists.
    int* current_out_degree_idx = (int*)calloc(N_VERTICES, sizeof(int));
    for (int i = 0; i < M_EDGES; i++) {
        int u = edge_list[i][0];
        int v = edge_list[i][1];
        adj[u][current_out_degree_idx[u]++] = v;
    }

    // Step 4: Free temporary structures.
    free(edge_list);
    free(current_out_degree_idx);
}

void run_computation() {
    // Kahn's algorithm for topological sorting (iterative).

    // The queue for vertices with in-degree 0.
    int* queue = (int*)malloc(N_VERTICES * sizeof(int));
    int head = 0, tail = 0;

    // We need to modify in-degrees, so we work on a copy.
    int* current_in_degree = (int*)malloc(N_VERTICES * sizeof(int));
    memcpy(current_in_degree, in_degree, N_VERTICES * sizeof(int));

    // Initialize the queue with all vertices having an in-degree of 0.
    for (int i = 0; i < N_VERTICES; i++) {
        if (current_in_degree[i] == 0) {
            queue[tail++] = i;
        }
    }

    int result_idx = 0;
    while (head < tail) {
        // Dequeue vertex u
        int u = queue[head++];
        topological_sort_result[result_idx++] = u;

        // For each neighbor v of u
        for (int i = 0; i < out_degree[u]; i++) {
            int v = adj[u][i];
            current_in_degree[v]--;

            // If in-degree of v becomes 0, enqueue it
            if (current_in_degree[v] == 0) {
                queue[tail++] = v;
            }
        }
    }

    // Free memory used within computation
    free(queue);
    free(current_in_degree);
    
    // Generate a final result to prevent dead code elimination
    final_result = 0;
    // By construction, the sort should succeed. No need to check result_idx.
    for (int i = 0; i < N_VERTICES; i++) {
        final_result += topological_sort_result[i];
    }
}

void cleanup() {
    for (int i = 0; i < N_VERTICES; i++) {
        free(adj[i]);
    }
    free(adj);
    free(out_degree);
    free(in_degree);
    free(topological_sort_result);
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
    printf("%lld\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
