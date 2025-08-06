#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>
#include <stdbool.h>

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

// --- BENCHMARK DATA AND PARAMETERS ---
#define NUM_LABELS 8 // Number of atomic propositions

typedef struct {
    int num_states_model;
    int num_states_automaton;
    float edge_density;
    uint32_t seed;

    // Model (Kripke Structure)
    int** model_adj;          // Adjacency list for model transitions
    int* model_adj_size;      // Number of successors for each model state
    int* model_labels;        // Label for each model state (atomic proposition)

    // Automaton (Büchi)
    int** automaton_trans;    // automaton_trans[state][label] -> next_state
    bool* is_accepting;       // Flag for each automaton state

    // Result
    int cycle_count;
} BenchmarkData;

static BenchmarkData g_data;

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_states_model> <num_states_automaton> <edge_density> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_states_model = atoi(argv[1]);
    g_data.num_states_automaton = atoi(argv[2]);
    g_data.edge_density = atof(argv[3]);
    g_data.seed = (uint32_t)atoi(argv[4]);
    g_data.cycle_count = 0;

    mt_seed(g_data.seed);

    // Allocate and generate model data
    g_data.model_adj = malloc(g_data.num_states_model * sizeof(int*));
    g_data.model_adj_size = malloc(g_data.num_states_model * sizeof(int));
    g_data.model_labels = malloc(g_data.num_states_model * sizeof(int));

    for (int i = 0; i < g_data.num_states_model; ++i) {
        g_data.model_labels[i] = mt_rand() % NUM_LABELS;

        int degree = (int)(g_data.num_states_model * g_data.edge_density);
        if (degree == 0 && g_data.edge_density > 0.0f) degree = 1;
        if (degree >= g_data.num_states_model) degree = g_data.num_states_model -1;

        g_data.model_adj_size[i] = degree;
        g_data.model_adj[i] = malloc(degree * sizeof(int));
        for (int j = 0; j < degree; ++j) {
            g_data.model_adj[i][j] = mt_rand() % g_data.num_states_model;
        }
    }

    // Allocate and generate Büchi automaton data
    g_data.is_accepting = malloc(g_data.num_states_automaton * sizeof(bool));
    g_data.automaton_trans = malloc(g_data.num_states_automaton * sizeof(int*));
    for (int i = 0; i < g_data.num_states_automaton; ++i) {
        g_data.is_accepting[i] = (mt_rand() % 4 == 0); // ~25% are accepting states
        g_data.automaton_trans[i] = malloc(NUM_LABELS * sizeof(int));
        for (int j = 0; j < NUM_LABELS; ++j) {
            g_data.automaton_trans[i][j] = mt_rand() % g_data.num_states_automaton;
        }
    }
}

void run_computation() {
    int V_prod = g_data.num_states_model * g_data.num_states_automaton;

    int* stack = malloc(V_prod * sizeof(int) * 2);
    int* accepting_states_to_check = malloc(V_prod * sizeof(int) * 2);
    uint8_t* visited = calloc(V_prod, sizeof(uint8_t));

    // Pass 1: Find all reachable accepting product states using DFS
    int stack_top = 0;
    int accepting_count = 0;
    int initial_prod_idx = 0; // (s=0, q=0)

    stack[stack_top++] = 0; // q
    stack[stack_top++] = 0; // s
    visited[initial_prod_idx] = 1;

    while (stack_top > 0) {
        int s = stack[--stack_top];
        int q = stack[--stack_top];

        if (g_data.is_accepting[q]) {
            accepting_states_to_check[accepting_count++] = q;
            accepting_states_to_check[accepting_count++] = s;
        }

        int label = g_data.model_labels[s];
        int q_next = g_data.automaton_trans[q][label];
        for (int i = 0; i < g_data.model_adj_size[s]; ++i) {
            int s_next = g_data.model_adj[s][i];
            int next_prod_idx = s_next * g_data.num_states_automaton + q_next;
            if (!visited[next_prod_idx]) {
                visited[next_prod_idx] = 1;
                stack[stack_top++] = q_next;
                stack[stack_top++] = s_next;
            }
        }
    }
    free(visited);

    // Pass 2: For each reachable accepting state, check if it's on a cycle
    uint8_t* inner_visited = malloc(V_prod * sizeof(uint8_t));

    for (int i = 0; i < accepting_count; i += 2) {
        int start_q = accepting_states_to_check[i];
        int start_s = accepting_states_to_check[i + 1];
        int target_prod_idx = start_s * g_data.num_states_automaton + start_q;

        memset(inner_visited, 0, V_prod * sizeof(uint8_t));
        stack_top = 0;
        stack[stack_top++] = start_q;
        stack[stack_top++] = start_s;
        inner_visited[target_prod_idx] = 1;

        bool cycle_found = false;
        while (stack_top > 0) {
            int s = stack[--stack_top];
            int q = stack[--stack_top];

            int label = g_data.model_labels[s];
            int q_next = g_data.automaton_trans[q][label];

            for (int j = 0; j < g_data.model_adj_size[s]; ++j) {
                int s_next = g_data.model_adj[s][j];
                int next_prod_idx = s_next * g_data.num_states_automaton + q_next;

                if (next_prod_idx == target_prod_idx) {
                    cycle_found = true;
                    break;
                }
                if (!inner_visited[next_prod_idx]) {
                    inner_visited[next_prod_idx] = 1;
                    stack[stack_top++] = q_next;
                    stack[stack_top++] = s_next;
                }
            }
            if (cycle_found) break;
        }

        if (cycle_found) {
            g_data.cycle_count++;
        }
    }

    free(stack);
    free(accepting_states_to_check);
    free(inner_visited);
}

void cleanup() {
    for (int i = 0; i < g_data.num_states_model; ++i) {
        free(g_data.model_adj[i]);
    }
    free(g_data.model_adj);
    free(g_data.model_adj_size);
    free(g_data.model_labels);

    for (int i = 0; i < g_data.num_states_automaton; ++i) {
        free(g_data.automaton_trans[i]);
    }
    free(g_data.automaton_trans);
    free(g_data.is_accepting);
}

int main(int argc, char* argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%d\n", g_data.cycle_count);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
