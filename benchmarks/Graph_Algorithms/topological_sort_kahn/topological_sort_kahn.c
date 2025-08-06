#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

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

// --- Global Data Structures and Parameters ---
int num_vertices;
int num_edges;

// Adjacency list node
typedef struct AdjListNode {
    int dest;
    struct AdjListNode* next;
} AdjListNode;

// Graph representation
AdjListNode** adj_lists;
int* in_degree;

// Data for computation and result
int* top_order;
int* queue;
long long final_result;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_vertices> <num_edges> <seed>\n", argv[0]);
        exit(1);
    }

    num_vertices = atoi(argv[1]);
    num_edges = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (num_vertices <= 1 || num_edges < 0) {
        fprintf(stderr, "Invalid arguments: num_vertices > 1 and num_edges >= 0\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate memory for graph data structures
    adj_lists = (AdjListNode**)malloc(num_vertices * sizeof(AdjListNode*));
    in_degree = (int*)calloc(num_vertices, sizeof(int));
    top_order = (int*)malloc(num_vertices * sizeof(int));
    queue = (int*)malloc(num_vertices * sizeof(int));

    if (!adj_lists || !in_degree || !top_order || !queue) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    for (int i = 0; i < num_vertices; i++) {
        adj_lists[i] = NULL;
    }

    // Generate a Directed Acyclic Graph (DAG)
    // To guarantee a DAG, we only add edges (u, v) where u < v.
    for (int i = 0; i < num_edges; i++) {
        uint32_t u, v;
        do {
            u = mt_rand() % num_vertices;
            v = mt_rand() % num_vertices;
        } while (u == v); // Avoid self-loops

        if (u > v) { // Ensure u < v
            uint32_t temp = u;
            u = v;
            v = temp;
        }

        // Add edge from u to v
        AdjListNode* newNode = (AdjListNode*)malloc(sizeof(AdjListNode));
        if (!newNode) {
            fprintf(stderr, "FATAL: Mid-generation memory allocation failed.\n");
            exit(1);
        }
        newNode->dest = v;
        newNode->next = adj_lists[u];
        adj_lists[u] = newNode;
        in_degree[v]++;
    }
}

void run_computation() {
    int queue_head = 0;
    int queue_tail = 0;

    // 1. Initialize the queue with all vertices having an in-degree of 0
    for (int i = 0; i < num_vertices; i++) {
        if (in_degree[i] == 0) {
            queue[queue_tail++] = i;
        }
    }

    int top_order_idx = 0;

    // 2. Process vertices from the queue (Kahn's algorithm)
    while (queue_head < queue_tail) {
        // Dequeue vertex u
        int u = queue[queue_head++];
        top_order[top_order_idx++] = u;

        // Iterate over all adjacent vertices of the dequeued vertex u
        AdjListNode* current_node = adj_lists[u];
        while (current_node != NULL) {
            int v = current_node->dest;
            // Reduce in-degree of v. If it becomes 0, add it to the queue.
            in_degree[v]--;
            if (in_degree[v] == 0) {
                queue[queue_tail++] = v;
            }
            current_node = current_node->next;
        }
    }

    // 3. Compute a checksum of the topological order to prevent dead code elimination.
    // Our graph generation guarantees a valid topological sort (no cycles).
    final_result = 0;
    for (int i = 0; i < top_order_idx; i++) {
        final_result += (long long)top_order[i] * (i + 1);
    }
}

void cleanup() {
    // Free adjacency list nodes
    for (int i = 0; i < num_vertices; i++) {
        AdjListNode* current = adj_lists[i];
        while (current != NULL) {
            AdjListNode* temp = current;
            current = current->next;
            free(temp);
        }
    }
    // Free main data structure arrays
    free(adj_lists);
    free(in_degree);
    free(top_order);
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

    // Print the final computational result to stdout
    printf("%lld\n", final_result);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
