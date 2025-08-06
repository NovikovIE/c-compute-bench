#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <limits.h>
#include <stdbool.h>

// --- START MERSENNE TWISTER --- 
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

// Benchmark parameters and data structures
int V; // Number of vertices
int E; // Number of edges

// Adjacency list representation of the graph
struct Edge {
    int to;
    int weight;
    struct Edge* next;
};
struct Edge** adj_list;

// Arrays for Prim's algorithm
int* key;      // Key values used to pick minimum weight edge
bool* in_mst;  // To represent set of vertices included in MST

// Final result
long long total_mst_weight;

// Helper to add an undirected edge
void add_edge(int u, int v, int weight) {
    struct Edge* edge1 = (struct Edge*)malloc(sizeof(struct Edge));
    edge1->to = v;
    edge1->weight = weight;
    edge1->next = adj_list[u];
    adj_list[u] = edge1;

    struct Edge* edge2 = (struct Edge*)malloc(sizeof(struct Edge));
    edge2->to = u;
    edge2->weight = weight;
    edge2->next = adj_list[v];
    adj_list[v] = edge2;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_vertices> <num_edges> <seed>\n", argv[0]);
        exit(1);
    }

    V = atoi(argv[1]);
    E = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (V <= 0 || E < V - 1) {
        fprintf(stderr, "Error: num_vertices must be > 0 and num_edges >= num_vertices - 1 for a connected graph.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate adjacency list and initialize to NULL
    adj_list = (struct Edge**)calloc(V, sizeof(struct Edge*));
    if (!adj_list) {
        perror("Failed to allocate memory for adjacency list");
        exit(1);
    }

    // To guarantee connectivity, first create a path graph 0-1-2-...-V-1
    for (int i = 0; i < V - 1; ++i) {
        int weight = (mt_rand() % 1000) + 1;
        add_edge(i, i + 1, weight);
    }

    // Add remaining edges randomly to create a denser graph
    int edges_to_add = E - (V - 1);
    for (int i = 0; i < edges_to_add; ++i) {
        int u, v;
        do {
            u = mt_rand() % V;
            v = mt_rand() % V;
        } while (u == v); // Avoid self-loops
        int weight = (mt_rand() % 1000) + 1;
        add_edge(u, v, weight);
    }

    // Allocate memory for Prim's algorithm arrays
    key = (int*)malloc(V * sizeof(int));
    in_mst = (bool*)malloc(V * sizeof(bool));
    if (!key || !in_mst) {
        perror("Failed to allocate memory for Prim's arrays");
        exit(1);
    }
}

// Helper to find the vertex with minimum key value, from the set of vertices
// not yet included in the MST. This is the O(V) scan part.
int min_key_vertex() {
    int min = INT_MAX, min_index = -1;
    for (int v = 0; v < V; v++) {
        if (!in_mst[v] && key[v] < min) {
            min = key[v];
            min_index = v;
        }
    }
    return min_index;
}

void run_computation() {
    // Initialize all keys as infinite and in_mst as false
    for (int i = 0; i < V; i++) {
        key[i] = INT_MAX;
        in_mst[i] = false;
    }

    // Start with the first vertex.
    key[0] = 0;

    // The MST will have V vertices. Loop V-1 times to find edges.
    for (int count = 0; count < V; count++) {
        // Pick the minimum key vertex from the set of vertices not yet in MST
        int u = min_key_vertex();

        // If u is -1, it means the graph is disconnected. Our setup prevents this.
        if (u == -1) break;

        // Add the picked vertex to the MST Set
        in_mst[u] = true;

        // Update key value of all adjacent vertices of the picked vertex.
        struct Edge* edge = adj_list[u];
        while (edge != NULL) {
            int v = edge->to;
            int weight = edge->weight;
            if (!in_mst[v] && weight < key[v]) {
                key[v] = weight;
            }
            edge = edge->next;
        }
    }

    // Calculate the total weight of the MST by summing up the keys
    total_mst_weight = 0;
    for (int i = 0; i < V; i++) {
        if (key[i] != INT_MAX) { // Should be all vertices in a connected graph
             total_mst_weight += key[i];
        }
    }
}

void cleanup() {
    for (int i = 0; i < V; i++) {
        struct Edge* current = adj_list[i];
        while (current != NULL) {
            struct Edge* next = current->next;
            free(current);
            current = next;
        }
    }
    free(adj_list);
    free(key);
    free(in_mst);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%lld\n", total_mst_weight);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
