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

// --- Benchmark Specific Code ---

// Parameters
int NUM_VERTICES;
long NUM_EDGES;

// Adjacency list node
typedef struct Node {
    int vertex;
    struct Node* next;
} Node;

// Global data structures
Node** adj_list; // Adjacency list representation of the graph
int* visited;    // For tracking visited vertices during traversal
int* stack;      // Explicit stack for iterative DFS to prevent stack overflow
int final_result; // Stores the final count of connected components

// Helper to add an edge to an undirected graph
void add_edge(int src, int dest) {
    // Add an edge from src to dest
    Node* newNodeSrc = (Node*)malloc(sizeof(Node));
    newNodeSrc->vertex = dest;
    newNodeSrc->next = adj_list[src];
    adj_list[src] = newNodeSrc;

    // Add an edge from dest to src
    Node* newNodeDest = (Node*)malloc(sizeof(Node));
    newNodeDest->vertex = src;
    newNodeDest->next = adj_list[dest];
    adj_list[dest] = newNodeDest;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_vertices> <num_edges> <seed>\n", argv[0]);
        exit(1);
    }

    NUM_VERTICES = atoi(argv[1]);
    NUM_EDGES = atol(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (NUM_VERTICES <= 0 || NUM_EDGES < 0) {
        fprintf(stderr, "Error: Number of vertices must be positive, and edges non-negative.\n");
        exit(1);
    }

    mt_seed(seed);

    adj_list = (Node**)calloc(NUM_VERTICES, sizeof(Node*));
    visited = (int*)malloc(NUM_VERTICES * sizeof(int));
    stack = (int*)malloc(NUM_VERTICES * sizeof(int));

    if (!adj_list || !visited || !stack) {
        fprintf(stderr, "Failed to allocate memory for graph structures.\n");
        exit(1);
    }

    // Generate random edges
    for (long i = 0; i < NUM_EDGES; i++) {
        int u, v;
        do {
            u = mt_rand() % NUM_VERTICES;
            v = mt_rand() % NUM_VERTICES;
        } while (u == v); // Avoid self-loops
        add_edge(u, v);
    }
}

void run_computation() {
    int component_count = 0;
    
    // Initialize visited array to 0 (false)
    memset(visited, 0, NUM_VERTICES * sizeof(int));

    // Iterate through all vertices
    for (int i = 0; i < NUM_VERTICES; ++i) {
        if (!visited[i]) {
            // Found a new unvisited vertex, so it's a new component
            component_count++;
            
            // Iterative DFS to traverse the component
            int stack_top = -1;
            stack[++stack_top] = i;
            visited[i] = 1;

            while (stack_top != -1) {
                int u = stack[stack_top--];
                
                // Visit all adjacent vertices
                Node* current = adj_list[u];
                while (current != NULL) {
                    int v = current->vertex;
                    if (!visited[v]) {
                        visited[v] = 1;
                        stack[++stack_top] = v;
                    }
                    current = current->next;
                }
            }
        }
    }
    final_result = component_count;
}

void cleanup() {
    for (int i = 0; i < NUM_VERTICES; i++) {
        Node* current = adj_list[i];
        while (current != NULL) {
            Node* temp = current;
            current = current->next;
            free(temp);
        }
    }
    free(adj_list);
    free(visited);
    free(stack);
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
    printf("%d\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
