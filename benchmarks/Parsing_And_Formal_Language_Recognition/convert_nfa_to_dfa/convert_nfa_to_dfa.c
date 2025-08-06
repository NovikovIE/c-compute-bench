#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

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

// ----------- Global Data Structures -----------

// Parameters
int P_NUM_NFA_STATES;
int P_ALPHABET_SIZE;

// NFA representation (input)
uint64_t **g_nfa_transitions; // [state][symbol] -> bitmask of next states
uint64_t g_nfa_accept_states; // bitmask of NFA accept states

// --- Data for run_computation and cleanup ---

// A map from an NFA state set (uint64_t) to a DFA state index (int).
typedef struct {
    uint64_t key;
    int value_plus_one; // Using 0 as empty, so stored value is index + 1
} MapEntry;

MapEntry *g_dfa_map;
size_t g_dfa_map_capacity;

// Worklist for subset construction algorithm (a simple queue)
int *g_worklist;
int g_worklist_head;
int g_worklist_tail;

// Consolidated DFA representation
uint64_t *g_dfa_state_sets;
int **g_dfa_transitions;
char *g_dfa_accept_states;
int g_dfa_states_count;
size_t g_dfa_capacity;

// Final result
int g_final_result;

// ----------- Helper Functions for Hash Map -----------

static void map_insert(uint64_t key, int dfa_idx) {
    size_t h = key % g_dfa_map_capacity;
    while (g_dfa_map[h].value_plus_one != 0) {
        h = (h + 1) % g_dfa_map_capacity;
    }
    g_dfa_map[h].key = key;
    g_dfa_map[h].value_plus_one = dfa_idx + 1;
}

static int map_find(uint64_t key) {
    size_t h = key % g_dfa_map_capacity;
    while (g_dfa_map[h].value_plus_one != 0) {
        if (g_dfa_map[h].key == key) {
            return g_dfa_map[h].value_plus_one - 1;
        }
        h = (h + 1) % g_dfa_map_capacity;
    }
    return -1; // Not found
}

// ----------- Benchmark Functions -----------

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s num_nfa_states alphabet_size seed\n", argv[0]);
        exit(1);
    }
    P_NUM_NFA_STATES = atoi(argv[1]);
    P_ALPHABET_SIZE = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    if (P_NUM_NFA_STATES <= 0 || P_NUM_NFA_STATES > 64) {
        fprintf(stderr, "FATAL: num_nfa_states must be between 1 and 64.\n");
        exit(1);
    }

    mt_seed(seed);

    // 1. Allocate and generate the NFA transition table
    g_nfa_transitions = (uint64_t **)malloc(P_NUM_NFA_STATES * sizeof(uint64_t *));
    if (!g_nfa_transitions) {
        perror("malloc failed");
        exit(1);
    }
    for (int i = 0; i < P_NUM_NFA_STATES; i++) {
        g_nfa_transitions[i] = (uint64_t *)malloc(P_ALPHABET_SIZE * sizeof(uint64_t));
        if (!g_nfa_transitions[i]) {
            perror("malloc failed");
            exit(1);
        }
        for (int j = 0; j < P_ALPHABET_SIZE; j++) {
            uint64_t transition_mask = 0;
            int num_transitions = 1 + (mt_rand() % 3); // 1 to 3 transitions on average
            for (int k = 0; k < num_transitions; k++) {
                transition_mask |= (1ULL << (mt_rand() % P_NUM_NFA_STATES));
            }
            g_nfa_transitions[i][j] = transition_mask;
        }
    }

    // 2. Generate NFA accepting states
    g_nfa_accept_states = 0;
    for (int i = 0; i < P_NUM_NFA_STATES; i++) {
        if (mt_rand() % 5 == 0) { // ~20% of states are accepting
            g_nfa_accept_states |= (1ULL << i);
        }
    }
    if (g_nfa_accept_states == 0 && P_NUM_NFA_STATES > 0) {
        g_nfa_accept_states |= (1ULL << (mt_rand() % P_NUM_NFA_STATES));
    }

    // 3. Pre-allocate all memory for the DFA construction
    g_dfa_capacity = 1ULL << P_NUM_NFA_STATES;
    
    g_dfa_state_sets = (uint64_t *)malloc(g_dfa_capacity * sizeof(uint64_t));
    g_dfa_accept_states = (char *)malloc(g_dfa_capacity * sizeof(char));
    g_worklist = (int *)malloc(g_dfa_capacity * sizeof(int));
    
    g_dfa_transitions = (int **)malloc(g_dfa_capacity * sizeof(int *));
    for (size_t i = 0; i < g_dfa_capacity; i++) {
        g_dfa_transitions[i] = (int *)malloc(P_ALPHABET_SIZE * sizeof(int));
    }

    // Set map capacity to 2x max states to keep load factor <= 0.5
    g_dfa_map_capacity = g_dfa_capacity * 2;
    g_dfa_map = (MapEntry *)calloc(g_dfa_map_capacity, sizeof(MapEntry));

    // Check for allocation failures
    if (!g_dfa_state_sets || !g_dfa_accept_states || !g_worklist || !g_dfa_transitions || !g_dfa_map) {
        fprintf(stderr, "FATAL: Failed to allocate memory for DFA construction.\n");
        exit(1);
    }
}

void run_computation() {
    g_worklist_head = 0;
    g_worklist_tail = 0;
    g_dfa_states_count = 0;

    // State 0 is the trap state (for empty NFA sets)
    uint64_t empty_set = 0;
    g_dfa_state_sets[g_dfa_states_count] = empty_set;
    g_dfa_accept_states[g_dfa_states_count] = 0;
    for (int c = 0; c < P_ALPHABET_SIZE; c++) {
        g_dfa_transitions[g_dfa_states_count][c] = 0; // Trap state transitions to itself
    }
    map_insert(empty_set, g_dfa_states_count);
    g_dfa_states_count++;

    // The initial DFA state is the set containing the NFA start state (state 0)
    uint64_t start_set = 1ULL << 0;
    int start_dfa_idx = map_find(start_set);
    if (start_dfa_idx == -1) {
        start_dfa_idx = g_dfa_states_count;
        g_dfa_state_sets[start_dfa_idx] = start_set;
        map_insert(start_set, start_dfa_idx);
        g_worklist[g_worklist_tail++] = start_dfa_idx;
        g_dfa_states_count++;
    }

    // Main subset construction loop
    while (g_worklist_head < g_worklist_tail) {
        int current_dfa_idx = g_worklist[g_worklist_head++];
        uint64_t current_nfa_set = g_dfa_state_sets[current_dfa_idx];

        g_dfa_accept_states[current_dfa_idx] = ((current_nfa_set & g_nfa_accept_states) != 0);
        
        for (int c = 0; c < P_ALPHABET_SIZE; c++) {
            uint64_t next_nfa_set = 0;
            for (int i = 0; i < P_NUM_NFA_STATES; i++) {
                if ((current_nfa_set >> i) & 1) {
                    next_nfa_set |= g_nfa_transitions[i][c];
                }
            }

            int next_dfa_idx = map_find(next_nfa_set);
            if (next_dfa_idx == -1) {
                next_dfa_idx = g_dfa_states_count;
                g_dfa_state_sets[next_dfa_idx] = next_nfa_set;
                map_insert(next_nfa_set, next_dfa_idx);
                g_worklist[g_worklist_tail++] = next_dfa_idx;
                g_dfa_states_count++;
            }
            g_dfa_transitions[current_dfa_idx][c] = next_dfa_idx;
        }
    }

    g_final_result = g_dfa_states_count;
}

void cleanup() {
    // Free NFA structures
    for (int i = 0; i < P_NUM_NFA_STATES; i++) {
        free(g_nfa_transitions[i]);
    }
    free(g_nfa_transitions);

    // Free DFA structures
    free(g_dfa_state_sets);
    free(g_dfa_accept_states);
    free(g_worklist);
    for (size_t i = 0; i < g_dfa_capacity; i++) {
        free(g_dfa_transitions[i]);
    }
    free(g_dfa_transitions);
    free(g_dfa_map);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%d\n", g_final_result);
    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}
