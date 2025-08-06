#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- Mersenne Twister (MT19937) psychic-disco-prng --- 
// --- DO NOT MODIFY THIS SECTION ----------------------
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
// --- END OF MT19937 SECTION ------------------------

// --- Benchmark Globals and Types ---
// Adjacency list node for graph representation
typedef struct AdjListNode {
    int dest;
    struct AdjListNode* next;
} AdjListNode;

// Graph parameters and data structures
int num_vertices;
int num_edges;
AdjListNode** adj_list; // Adjacency list representation of the graph

// State for the computation
enum { WHITE, GRAY, BLACK };
int* color; // Color array for DFS
volatile int final_result; // Accumulated result to prevent dead-code elimination

// --- Cycle Detection Helper --- 
// Recursive utility function for DFS-based cycle detection.
// It counts back edges, which indicate cycles.
static void detect_cycle_util(int u) {
    color[u] = GRAY; // Mark current node as visiting (in recursion stack)

    AdjListNode* node = adj_list[u];
    while (node != NULL) {
        int v = node->dest;

        // If an adjacent vertex 'v' is gray, then there is a back edge.
        if (color[v] == GRAY) {
            final_result++;
        }

        // If 'v' has not been visited yet, recurse on it.
        if (color[v] == WHITE) {
            detect_cycle_util(v);
        }

        node = node->next;
    }

    // Mark current node as visited (finished)
    color[u] = BLACK;
}

// --- Benchmark Functions ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_vertices> <num_edges> <seed>\n", argv[0]);
        exit(1);
    }

    num_vertices = atoi(argv[1]);
    num_edges = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (num_vertices <= 0 || num_edges < 0) {
        fprintf(stderr, "FATAL: Number of vertices and edges must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate adjacency list
    adj_list = (AdjListNode**)malloc(num_vertices * sizeof(AdjListNode*));
    if (!adj_list) {
        fprintf(stderr, "FATAL: Memory allocation failed for adjacency list.\n");
        exit(1);
    }
    for (int i = 0; i < num_vertices; i++) {
        adj_list[i] = NULL;
    }

    // Generate random edges (u -> v)
    for (int i = 0; i < num_edges; i++) {
        int u = mt_rand() % num_vertices;
        int v = mt_rand() % num_vertices;

        // Create and add new node to the front of the adjacency list for u
        AdjListNode* newNode = (AdjListNode*)malloc(sizeof(AdjListNode));
        if (!newNode) {
            fprintf(stderr, "FATAL: Memory allocation failed for graph node.\n");
            exit(1);
        }
        newNode->dest = v;
        newNode->next = adj_list[u];
        adj_list[u] = newNode;
    }

    // Allocate color array for DFS
    color = (int*)malloc(num_vertices * sizeof(int));
    if (!color) {
        fprintf(stderr, "FATAL: Memory allocation failed for color array.\n");
        exit(1);
    }
}

void run_computation() {
    final_result = 0;

    // Initialize all vertices to WHITE (not visited)
    for (int i = 0; i < num_vertices; i++) {
        color[i] = WHITE;
    }

    // Perform DFS from each unvisited vertex to find cycles
    for (int i = 0; i < num_vertices; i++) {
        if (color[i] == WHITE) {
            detect_cycle_util(i);
        }
    }
}

void cleanup() {
    for (int i = 0; i < num_vertices; i++) {
        AdjListNode* current = adj_list[i];
        while (current != NULL) {
            AdjListNode* next = current->next;
            free(current);
            current = next;
        }
    }
    free(adj_list);
    free(color);
}

// --- Main Function ---
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%d\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
