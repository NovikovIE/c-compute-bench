#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

// Mersenne Twister (MT19937) PRNG
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

// Benchmark parameters
int num_vertices;
long long num_edges;

// Graph data structures
int** adj_list;      // Adjacency list for the graph
int* adj_list_len;   // Number of outgoing edges for each vertex
int* in_degree;      // Master copy of in-degrees for each vertex

// Data structures for computation
int* sorted_order;       // To store the result of the topological sort
int* queue;              // Queue for vertices with in-degree 0
int* in_degree_comp;     // Working copy of in-degrees for computation

// Final result
long long final_checksum = 0;

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_vertices> <num_edges> <seed>\n", argv[0]);
        exit(1);
    }
    num_vertices = atoi(argv[1]);
    num_edges = atoll(argv[2]);
    uint32_t seed = atoi(argv[3]);

    mt_seed(seed);

    // Allocate master data structures
    adj_list = (int**)malloc(num_vertices * sizeof(int*));
    adj_list_len = (int*)calloc(num_vertices, sizeof(int));
    in_degree = (int*)calloc(num_vertices, sizeof(int));
    if (!adj_list || !adj_list_len || !in_degree) {
        fprintf(stderr, "Failed to allocate graph metadata\n");
        exit(1);
    }

    // --- Graph Generation (produces a DAG) ---
    // Pass 1: Generate edges and count degrees to determine allocation sizes.
    int* temp_u = (int*)malloc(num_edges * sizeof(int));
    int* temp_v = (int*)malloc(num_edges * sizeof(int));
    if (!temp_u || !temp_v) {
        fprintf(stderr, "Failed to allocate temporary edge list\n");
        exit(1);
    }

    for (long long i = 0; i < num_edges; ++i) {
        int u, v;
        // Ensure u < v to create a Directed Acyclic Graph (DAG).
        do {
            u = mt_rand() % num_vertices;
            v = mt_rand() % num_vertices;
        } while (u >= v);
        
        temp_u[i] = u;
        temp_v[i] = v;
        adj_list_len[u]++;
        in_degree[v]++;
    }

    // Pass 2: Allocate the adjacency lists with the exact required sizes.
    for (int i = 0; i < num_vertices; ++i) {
        if (adj_list_len[i] > 0) {
            adj_list[i] = (int*)malloc(adj_list_len[i] * sizeof(int));
        } else {
            adj_list[i] = NULL;
        }
    }

    // Pass 3: Populate the adjacency lists.
    int* current_adj_idx = (int*)calloc(num_vertices, sizeof(int));
    for (long long i = 0; i < num_edges; ++i) {
        int u = temp_u[i];
        int v = temp_v[i];
        adj_list[u][current_adj_idx[u]++] = v;
    }

    free(temp_u);
    free(temp_v);
    free(current_adj_idx);

    // Allocate data structures needed for the computation phase
    sorted_order = (int*)malloc(num_vertices * sizeof(int));
    queue = (int*)malloc(num_vertices * sizeof(int));
    in_degree_comp = (int*)malloc(num_vertices * sizeof(int));
    if (!sorted_order || !queue || !in_degree_comp) {
        fprintf(stderr, "Failed to allocate computation structures\n");
        exit(1);
    }
}

void run_computation() {
    // Kahn's algorithm for topological sort
    memcpy(in_degree_comp, in_degree, num_vertices * sizeof(int));

    int queue_head = 0;
    int queue_tail = 0;

    // Initialize the queue with all vertices having an in-degree of 0
    for (int i = 0; i < num_vertices; ++i) {
        if (in_degree_comp[i] == 0) {
            queue[queue_tail++] = i;
        }
    }

    int sorted_count = 0;

    while (queue_head < queue_tail) {
        int u = queue[queue_head++];
        sorted_order[sorted_count++] = u;

        // For each neighbor of u, decrement its in-degree
        for (int i = 0; i < adj_list_len[u]; ++i) {
            int v = adj_list[u][i];
            in_degree_comp[v]--;
            if (in_degree_comp[v] == 0) {
                queue[queue_tail++] = v;
            }
        }
    }

    // This benchmark assumes a DAG, a check for cycles is omitted.
    // if (sorted_count < num_vertices) { // Graph has a cycle }

    // Calculate a checksum to prevent dead code elimination.
    long long checksum = 0;
    for (int i = 0; i < sorted_count; ++i) {
        checksum = (checksum * 31) + sorted_order[i];
    }
    final_checksum = checksum;
}

void cleanup() {
    for (int i = 0; i < num_vertices; ++i) {
        free(adj_list[i]);
    }
    free(adj_list);
    free(adj_list_len);
    free(in_degree);
    free(sorted_order);
    free(queue);
    free(in_degree_comp);
}

int main(int argc, char* argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout
    printf("%lld\n", final_checksum);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
