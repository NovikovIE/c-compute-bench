#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

// --- BEGIN MERSENNE TWISTER (DO NOT MODIFY) ---
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

// Represents a state in the Kripke structure (state-transition system)
typedef struct {
    int* transitions;
    // In LTL model checking, we look for cycles through states that violate a property.
    // We simplify this by flagging some states as "bad".
    int is_bad_state;
} State;

// Global struct to hold all benchmark data
struct {
    int num_states;
    int num_transitions_per_state;
    int ltl_formula_complexity;

    State* graph;
    int* visited_outer;      // For the main reachability traversal
    int* visited_inner;      // For the cycle detection sub-problem
    int* on_stack;           // For the cycle detection sub-problem
    int* checked_for_cycle;  // To prune redundant sub-problem searches
    
    // The result is the number of counterexamples (bad cycles) found.
    volatile int final_result;
} g_data;

// Forward declarations for recursive functions
void main_dfs(int u);
int cycle_check_dfs(int u);

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_states num_transitions_per_state ltl_formula_complexity seed\n", argv[0]);
        exit(1);
    }

    g_data.num_states = atoi(argv[1]);
    g_data.num_transitions_per_state = atoi(argv[2]);
    g_data.ltl_formula_complexity = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    mt_seed(seed);

    g_data.graph = (State*)malloc(g_data.num_states * sizeof(State));
    g_data.visited_outer = (int*)malloc(g_data.num_states * sizeof(int));
    g_data.visited_inner = (int*)malloc(g_data.num_states * sizeof(int));
    g_data.on_stack = (int*)malloc(g_data.num_states * sizeof(int));
    g_data.checked_for_cycle = (int*)malloc(g_data.num_states * sizeof(int));

    if (!g_data.graph || !g_data.visited_outer || !g_data.visited_inner || !g_data.on_stack || !g_data.checked_for_cycle) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    for (int i = 0; i < g_data.num_states; ++i) {
        g_data.graph[i].transitions = (int*)malloc(g_data.num_transitions_per_state * sizeof(int));
        if (!g_data.graph[i].transitions) {
            fprintf(stderr, "FATAL: Memory allocation failed for transitions.\n");
            exit(1);
        }

        for (int j = 0; j < g_data.num_transitions_per_state; ++j) {
            g_data.graph[i].transitions[j] = mt_rand() % g_data.num_states;
        }

        // The LTL formula complexity determines the probability of a state being "bad".
        // A property like F(G(!p)) requires finding a cycle of !p states.
        g_data.graph[i].is_bad_state = (mt_rand() % 100) < g_data.ltl_formula_complexity;
    }
}

// Inner DFS: checks for a cycle of bad states starting from u.
int cycle_check_dfs(int u) {
    g_data.visited_inner[u] = 1;
    g_data.on_stack[u] = 1;

    for (int i = 0; i < g_data.num_transitions_per_state; i++) {
        int v = g_data.graph[u].transitions[i];
        if (g_data.graph[v].is_bad_state) { // Only traverse through bad states
            if (!g_data.visited_inner[v]) {
                if (cycle_check_dfs(v)) return 1;
            } else if (g_data.on_stack[v]) {
                return 1; // Cycle detected
            }
        }
    }

    g_data.on_stack[u] = 0;
    return 0;
}

// Outer DFS: Traverses the graph to find reachable states.
void main_dfs(int u) {
    g_data.visited_outer[u] = 1;

    // If we find a bad state that we haven't already processed, 
    // start a sub-problem to check for a cycle from it.
    if (g_data.graph[u].is_bad_state && !g_data.checked_for_cycle[u]) {
        memset(g_data.visited_inner, 0, g_data.num_states * sizeof(int));
        memset(g_data.on_stack, 0, g_data.num_states * sizeof(int));
        
        if (cycle_check_dfs(u)) {
            g_data.final_result++;
        }

        // Mark all states visited in the inner search as checked to avoid redundant work.
        for (int i = 0; i < g_data.num_states; i++) {
            if (g_data.visited_inner[i]) {
                g_data.checked_for_cycle[i] = 1;
            }
        }
    }

    // Continue the main traversal.
    for (int i = 0; i < g_data.num_transitions_per_state; ++i) {
        int v = g_data.graph[u].transitions[i];
        if (!g_data.visited_outer[v]) {
            main_dfs(v);
        }
    }
}

// Simulates an explicit-state model checker trying to find a counterexample
// to an LTL formula. The core task is a nested depth-first search.
void run_computation() {
    g_data.final_result = 0;
    memset(g_data.visited_outer, 0, g_data.num_states * sizeof(int));
    memset(g_data.checked_for_cycle, 0, g_data.num_states * sizeof(int));

    // In a real model checker, we'd start from all initial states. 
    // For this benchmark, we start from state 0.
    main_dfs(0);
    
    // Ensure we check all disconnected components as well.
    for(int i = 0; i < g_data.num_states; ++i) {
        if (!g_data.visited_outer[i]) {
            main_dfs(i);
        }
    }
}

void cleanup() {
    for (int i = 0; i < g_data.num_states; ++i) {
        free(g_data.graph[i].transitions);
    }
    free(g_data.graph);
    free(g_data.visited_outer);
    free(g_data.visited_inner);
    free(g_data.on_stack);
    free(g_data.checked_for_cycle);
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
    printf("%d\n", g_data.final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
