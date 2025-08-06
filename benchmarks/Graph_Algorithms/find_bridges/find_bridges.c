#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>

#define MIN(a, b) (((a) < (b)) ? (a) : (b))

// --- Mersenne Twister (DO NOT MODIFY) ---
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
struct AdjListNode {
    int dest;
    struct AdjListNode* next;
};

// Stack frame for iterative DFS
struct StackFrame {
    int u;
    struct AdjListNode* edge_iter;
};

static int V, E;
static struct AdjListNode** adj;       // Adjacency list array
static struct AdjListNode* edge_pool; // Memory pool for all nodes

// Arrays for bridge finding algorithm
static int* disc;
static int* low;
static int* parent;
static bool* visited;

// The final result: number of bridges found
static int bridge_count;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_vertices> <num_edges> <seed>\n", argv[0]);
        exit(1);
    }

    V = atoi(argv[1]);
    E = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    mt_seed(seed);

    // Allocate memory for graph and algorithm arrays
    adj = (struct AdjListNode**)malloc(V * sizeof(struct AdjListNode*));
    edge_pool = (struct AdjListNode*)malloc(E * 2 * sizeof(struct AdjListNode));
    disc = (int*)malloc(V * sizeof(int));
    low = (int*)malloc(V * sizeof(int));
    parent = (int*)malloc(V * sizeof(int));
    visited = (bool*)malloc(V * sizeof(bool));

    if (!adj || !edge_pool || !disc || !low || !parent || !visited) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize adjacency lists to empty
    for (int i = 0; i < V; i++) {
        adj[i] = NULL;
    }

    // Generate graph with E edges
    int edge_pool_idx = 0;
    for (int i = 0; i < E; i++) {
        int u = mt_rand() % V;
        int v = mt_rand() % V;

        if (u == v) {
            i--; // try again
            continue;
        }

        // Add edge from u to v
        struct AdjListNode* u_node = &edge_pool[edge_pool_idx++];
        u_node->dest = v;
        u_node->next = adj[u];
        adj[u] = u_node;

        // Add edge from v to u
        struct AdjListNode* v_node = &edge_pool[edge_pool_idx++];
        v_node->dest = u;
        v_node->next = adj[v];
        adj[v] = v_node;
    }
}

void run_computation() {
    // Initialize arrays and counters for the algorithm
    for (int i = 0; i < V; i++) {
        parent[i] = -1;
        visited[i] = false;
    }   
    bridge_count = 0;
    int time = 0;

    // Stack for iterative DFS, allocated on the heap to avoid stack overflow
    struct StackFrame* stack = (struct StackFrame*)malloc(V * sizeof(struct StackFrame));
    if (stack == NULL) {
        fprintf(stderr, "FATAL: Failed to allocate stack for DFS.\n");
        exit(1);
    }

    // Find bridges in all DFS trees
    for (int i = 0; i < V; i++) {
        if (!visited[i]) {
            int stack_top = -1;

            // Start DFS from vertex i
            stack_top++;
            stack[stack_top].u = i;
            stack[stack_top].edge_iter = adj[i];

            visited[i] = true;
            disc[i] = low[i] = ++time;
            // parent[i] is already -1

            while (stack_top >= 0) {
                int u = stack[stack_top].u;
                struct AdjListNode* node = stack[stack_top].edge_iter;

                if (node != NULL) {
                    int v = node->dest;
                    // Advance iterator for the current frame before exploring child
                    stack[stack_top].edge_iter = node->next;

                    if (visited[v]) {
                        if (v != parent[u]) {
                            low[u] = MIN(low[u], disc[v]);
                        }
                    } else {
                        // Tree edge: "descend" into v
                        visited[v] = true;
                        parent[v] = u;
                        disc[v] = low[v] = ++time;

                        // Push v onto the stack
                        stack_top++;
                        stack[stack_top].u = v;
                        stack[stack_top].edge_iter = adj[v];
                    }
                } else {
                    // Finished exploring all edges from u, "return" from recursion
                    // The u we are processing is stack[stack_top].u
                    
                    // Update parent's low-link value, if parent exists
                    if (parent[u] != -1) {
                        int p = parent[u];
                        low[p] = MIN(low[p], low[u]);
                        if (low[u] > disc[p]) {
                            bridge_count++;
                        }
                    }
                    
                    // Pop from stack
                    stack_top--;
                }
            }
        }
    }
    
    free(stack);
}

void cleanup() {
    free(adj);
    free(edge_pool);
    free(disc);
    free(low);
    free(parent);
    free(visited);
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
    printf("%d\n", bridge_count);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
