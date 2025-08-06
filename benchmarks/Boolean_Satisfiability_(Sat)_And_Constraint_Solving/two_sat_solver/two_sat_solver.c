#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <stdbool.h>
#include <string.h> 

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
// --- End Mersenne Twister ---

// --- Benchmark Globals ---
int num_variables;
int num_clauses;
int num_nodes;
size_t num_edges;

// Implication graph G and its reverse G_rev
// Stored as flattened adjacency lists (CSR format)
int* adj;         // Stores edges for G
int* adj_rev;     // Stores edges for G_rev
size_t* adj_pos;     // Stores starting position of each node's edges in adj
size_t* adj_pos_rev; // Stores starting position for adj_rev

int final_result; // 1 if satisfiable, 0 otherwise

// --- Helper Functions ---

// Maps a literal to a graph node index.
// Variable v (1-based) becomes node 2*(v-1).
// Negation ~v becomes node 2*(v-1) + 1.
static inline int literal_to_node(int literal) {
    if (literal > 0) {
        return 2 * (literal - 1);
    } else {
        return 2 * (-literal - 1) + 1;
    }
}

// --- Core Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_variables> <num_clauses> <seed>\n", argv[0]);
        exit(1);
    }
    num_variables = atoi(argv[1]);
    num_clauses = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);
    
    mt_seed(seed);
    
    num_nodes = 2 * num_variables;
    num_edges = 2 * (size_t)num_clauses;
    
    // Step 1: Generate all clause implications and store temporarily
    // Each clause (a or b) -> (~a => b) and (~b => a)
    // We store the two directed edges {u, v} for each clause.
    int (*implications)[4] = malloc((size_t)num_clauses * sizeof(int[4]));
    if (!implications) {
        perror("Failed to allocate implications buffer");
        exit(1);
    }

    for (int i = 0; i < num_clauses; ++i) {
        int v1_idx = mt_rand() % num_variables;
        int v2_idx;
        do {
            v2_idx = mt_rand() % num_variables;
        } while (v1_idx == v2_idx);

        // Literals are 1-based, positive for var, negative for negation
        int lit1 = (v1_idx + 1) * (mt_rand() % 2 ? 1 : -1);
        int lit2 = (v2_idx + 1) * (mt_rand() % 2 ? 1 : -1);

        // Implication ~lit1 => lit2
        implications[i][0] = literal_to_node(-lit1);
        implications[i][1] = literal_to_node(lit2);
        // Implication ~lit2 => lit1
        implications[i][2] = literal_to_node(-lit2);
        implications[i][3] = literal_to_node(lit1);
    }

    // Step 2 & 3: Count degrees and create CSR index pointers (prefix sum)
    adj_pos = calloc(num_nodes + 1, sizeof(size_t));
    adj_pos_rev = calloc(num_nodes + 1, sizeof(size_t));
    if (!adj_pos || !adj_pos_rev) {
        perror("Failed to allocate CSR positions");
        exit(1);
    }

    for (int i = 0; i < num_clauses; ++i) {
        adj_pos[implications[i][0] + 1]++;
        adj_pos_rev[implications[i][1] + 1]++;
        adj_pos[implications[i][2] + 1]++;
        adj_pos_rev[implications[i][3] + 1]++;
    }

    for (int i = 0; i < num_nodes; ++i) {
        adj_pos[i+1] += adj_pos[i];
        adj_pos_rev[i+1] += adj_pos_rev[i];
    }
    
    // Step 4: Allocate adjacency lists and fill them
    adj = malloc(num_edges * sizeof(int));
    adj_rev = malloc(num_edges * sizeof(int));
    if (!adj || !adj_rev) {
        perror("Failed to allocate adjacency lists");
        exit(1);
    }

    size_t* current_pos = malloc((num_nodes + 1) * sizeof(size_t));
    size_t* current_pos_rev = malloc((num_nodes + 1) * sizeof(size_t));
    if (!current_pos || !current_pos_rev) {
        perror("Failed to allocate CSR temp counters");
        exit(1);
    }
    memcpy(current_pos, adj_pos, (num_nodes + 1) * sizeof(size_t));
    memcpy(current_pos_rev, adj_pos_rev, (num_nodes + 1) * sizeof(size_t));

    for (int i = 0; i < num_clauses; ++i) {
        int u1 = implications[i][0];
        int v1 = implications[i][1];
        adj[current_pos[u1]++] = v1;
        adj_rev[current_pos_rev[v1]++] = u1;

        int u2 = implications[i][2];
        int v2 = implications[i][3];
        adj[current_pos[u2]++] = v2;
        adj_rev[current_pos_rev[v2]++] = u2;
    }

    free(implications);
    free(current_pos);
    free(current_pos_rev);
}


// --- Kosaraju's Algorithm components (used in run_computation) ---
static bool* visited;
static int* finish_stack;
static int finish_stack_top;
static int* scc_id;

// Iterative DFS for post-order traversal (part 1 of Kosaraju's)
// Simulates recursion to avoid stack overflow on large graphs.
static void iterative_dfs1(int start_node, int* dfs_stack, size_t* edge_idx_stack) {
    int dfs_stack_top = 0;

    // Push start node onto the traversal stack
    dfs_stack[dfs_stack_top] = start_node;
    edge_idx_stack[dfs_stack_top] = adj_pos[start_node];
    dfs_stack_top++;
    visited[start_node] = true;

    while (dfs_stack_top > 0) {
        int u = dfs_stack[dfs_stack_top - 1];
        size_t* current_edge_pos = &edge_idx_stack[dfs_stack_top - 1];

        bool found_new_node = false;
        while (*current_edge_pos < adj_pos[u + 1]) {
            int v = adj[*current_edge_pos];
            (*current_edge_pos)++;
            if (!visited[v]) {
                visited[v] = true;
                // Push v onto the stack
                dfs_stack[dfs_stack_top] = v;
                edge_idx_stack[dfs_stack_top] = adj_pos[v];
                dfs_stack_top++;
                found_new_node = true;
                break;
            }
        }

        // If all neighbors visited, u is finished. Pop it and add to finish_stack.
        if (!found_new_node) {
            finish_stack[finish_stack_top++] = u;
            dfs_stack_top--;
        }
    }
}

// Iterative DFS for pre-order traversal (part 2 of Kosaraju's)
static void iterative_dfs2(int start_node, int current_scc, int* dfs_stack) {
    int dfs_stack_top = 0;

    // Push start node
    dfs_stack[dfs_stack_top++] = start_node;
    visited[start_node] = true;
    scc_id[start_node] = current_scc;

    while (dfs_stack_top > 0) {
        int u = dfs_stack[--dfs_stack_top];

        for (size_t i = adj_pos_rev[u]; i < adj_pos_rev[u + 1]; ++i) {
            int v = adj_rev[i];
            if (!visited[v]) {
                visited[v] = true;
                scc_id[v] = current_scc;
                dfs_stack[dfs_stack_top++] = v;
            }
        }
    }
}

void run_computation() {
    // Allocate temporary arrays for computation
    visited = malloc(num_nodes * sizeof(bool));
    finish_stack = malloc(num_nodes * sizeof(int));
    scc_id = malloc(num_nodes * sizeof(int));

    // Allocate stacks for iterative DFS to prevent stack overflow
    int* dfs_stack = malloc(num_nodes * sizeof(int));
    size_t* edge_idx_stack = malloc(num_nodes * sizeof(size_t));

    if (!visited || !finish_stack || !scc_id || !dfs_stack || !edge_idx_stack) {
        perror("Failed to allocate computation buffers");
        exit(1);
    }

    // Step 1 of Kosaraju's algorithm: DFS on G to get finishing times
    finish_stack_top = 0;
    memset(visited, 0, num_nodes * sizeof(bool));
    for (int i = 0; i < num_nodes; ++i) {
        if (!visited[i]) {
            iterative_dfs1(i, dfs_stack, edge_idx_stack);
        }
    }

    // Step 2 of Kosaraju's algorithm: DFS on G_rev to find SCCs
    memset(visited, 0, num_nodes * sizeof(bool));
    int current_scc = 0;
    while (finish_stack_top > 0) {
        int u = finish_stack[--finish_stack_top];
        if (!visited[u]) {
            iterative_dfs2(u, current_scc++, dfs_stack);
        }
    }
    
    // Step 3: Check for satisfiability
    final_result = 1; // Assume satisfiable
    for (int i = 0; i < num_variables; ++i) {
        // Variable i corresponds to node 2*i
        // Negation ~i corresponds to node 2*i + 1
        if (scc_id[2 * i] == scc_id[2 * i + 1]) {
            final_result = 0; // Unsatisfiable
            break;
        }
    }

    // Free temporary arrays
    free(visited);
    free(finish_stack);
    free(scc_id);
    free(dfs_stack);
    free(edge_idx_stack);
}

void cleanup() {
    free(adj);
    free(adj_rev);
    free(adj_pos);
    free(adj_pos_rev);
}


int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%d\n", final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
