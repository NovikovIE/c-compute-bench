#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <stdbool.h>

// --- Mersenne Twister (MT19937) --- DO NOT MODIFY ---
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
// --- End of Mersenne Twister ---

// --- Benchmark Globals ---

// Adjacency list node
struct Node {
    int dest;
    struct Node* next;
};

int num_vertices;
int num_edges;

// Graph G and its transpose G_T
struct Node** adj;
struct Node** adj_T;

// Auxiliary data structures for Kosaraju's algorithm
bool* visited;
int* finish_stack;
int stack_top;

// Final result: number of Strongly Connected Components
int final_result;

// --- Helper functions for the algorithm ---

// First DFS: fills the finish_stack in post-order traversal
void dfs1(int u) {
    visited[u] = true;
    struct Node* node = adj[u];
    while (node != NULL) {
        if (!visited[node->dest]) {
            dfs1(node->dest);
        }
        node = node->next;
    }
    finish_stack[stack_top++] = u;
}

// Second DFS: traverses the transpose graph to find SCCs
void dfs2(int u) {
    visited[u] = true;
    struct Node* node = adj_T[u];
    while (node != NULL) {
        if (!visited[node->dest]) {
            dfs2(node->dest);
        }
        node = node->next;
    }
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s num_vertices num_edges seed\n", argv[0]);
        exit(1);
    }

    num_vertices = atoi(argv[1]);
    num_edges = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (num_vertices <= 0 || num_edges < 0) {
        fprintf(stderr, "FATAL: num_vertices must be positive and num_edges non-negative.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate memory for graph and auxiliary structures
    adj = (struct Node**)calloc(num_vertices, sizeof(struct Node*));
    adj_T = (struct Node**)calloc(num_vertices, sizeof(struct Node*));
    visited = (bool*)malloc(num_vertices * sizeof(bool));
    finish_stack = (int*)malloc(num_vertices * sizeof(int));

    if (!adj || !adj_T || !visited || !finish_stack) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Generate random graph
    for (int i = 0; i < num_edges; ++i) {
        int u = mt_rand() % num_vertices;
        int v = mt_rand() % num_vertices;

        // Add edge u -> v to the original graph
        struct Node* newNode = (struct Node*)malloc(sizeof(struct Node));
        if (!newNode) { fprintf(stderr, "FATAL: Memory allocation failed.\n"); exit(1); }
        newNode->dest = v;
        newNode->next = adj[u];
        adj[u] = newNode;

        // Add edge v -> u to the transpose graph
        struct Node* newTNode = (struct Node*)malloc(sizeof(struct Node));
        if (!newTNode) { fprintf(stderr, "FATAL: Memory allocation failed.\n"); exit(1); }
        newTNode->dest = u;
        newTNode->next = adj_T[v];
        adj_T[v] = newTNode;
    }
}

void run_computation() {
    // Step 1: Perform first DFS on G to get finish times (post-order)
    for (int i = 0; i < num_vertices; ++i) {
        visited[i] = false;
    }
    stack_top = 0;
    for (int i = 0; i < num_vertices; ++i) {
        if (!visited[i]) {
            dfs1(i);
        }
    }

    // Step 2: Perform second DFS on G_T in the order of finish times
    for (int i = 0; i < num_vertices; ++i) {
        visited[i] = false;
    }
    int scc_count = 0;
    while (stack_top > 0) {
        int u = finish_stack[--stack_top];
        if (!visited[u]) {
            dfs2(u);
            scc_count++;
        }
    }

    final_result = scc_count;
}

void cleanup() {
    for (int i = 0; i < num_vertices; ++i) {
        struct Node* current = adj[i];
        while (current != NULL) {
            struct Node* temp = current;
            current = current->next;
            free(temp);
        }
        current = adj_T[i];
        while (current != NULL) {
            struct Node* temp = current;
            current = current->next;
            free(temp);
        }
    }
    free(adj);
    free(adj_T);
    free(visited);
    free(finish_stack);
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
    printf("%d\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
