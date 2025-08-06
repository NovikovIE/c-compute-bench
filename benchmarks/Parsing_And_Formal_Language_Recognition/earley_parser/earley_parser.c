#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// MT19937 Generator (verbatim)
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
// End MT19937

#define MAX_RULE_BODY 2
// Note: Symbols are integers.
// > 0: Terminal
// < 0: Non-terminal
// 0: Not used
#define START_SYMBOL -1

// Data structures for Earley Parser
typedef struct {
    int head; // Non-terminal LHS
    int body[MAX_RULE_BODY];
    int body_len;
} Rule;

typedef struct {
    int rule_idx;   // Index into the grammar->rules array
    int dot_pos;    // Position of the dot in the rule body
    int origin_pos; // Start position in the input string (j in [A -> a.B, j])
} State;

typedef struct {
    State* states;
    int count;
    int capacity;
} StateSet;

// Benchmark parameters
static int INPUT_TOKEN_COUNT;
static int GRAMMAR_RULE_COUNT;
static int NUM_NON_TERMINALS;
static int NUM_TERMINALS;

// Global data structures
static Rule* grammar = NULL;
static int* input_tokens = NULL;
static StateSet* chart = NULL;
static long long final_result = 0;

// Helper to add a state, preventing duplicates
void add_state_to_set(State state, StateSet* set) {
    for (int i = 0; i < set->count; i++) {
        if (set->states[i].rule_idx == state.rule_idx &&
            set->states[i].dot_pos == state.dot_pos &&
            set->states[i].origin_pos == state.origin_pos) {
            return; // State already exists, do not add.
        }
    }

    if (set->count >= set->capacity) {
        set->capacity = set->capacity == 0 ? 32 : set->capacity * 2;
        set->states = (State*)realloc(set->states, set->capacity * sizeof(State));
        if (!set->states) {
            perror("Failed to reallocate state set");
            exit(EXIT_FAILURE);
        }
    }
    set->states[set->count++] = state;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <input_token_count> <grammar_rule_count> <num_non_terminals> <num_terminals> <seed>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    INPUT_TOKEN_COUNT = atoi(argv[1]);
    GRAMMAR_RULE_COUNT = atoi(argv[2]);
    NUM_NON_TERMINALS = atoi(argv[3]);
    NUM_TERMINALS = atoi(argv[4]);
    uint32_t seed = atoi(argv[5]);

    mt_seed(seed);

    // --- Generate Grammar ---
    // The grammar will contain GRAMMAR_RULE_COUNT + 1 rules.
    // Rule 0 is reserved for the initial state: S' -> S (where S is START_SYMBOL)
    grammar = (Rule*)malloc((GRAMMAR_RULE_COUNT + 1) * sizeof(Rule));
    if (!grammar) { perror("malloc grammar"); exit(EXIT_FAILURE); }

    // Rule 0: S' -> S
    grammar[0].head = 0; // Dummy start symbol S'
    grammar[0].body[0] = START_SYMBOL;
    grammar[0].body_len = 1;

    // Generate other rules
    for (int i = 1; i <= GRAMMAR_RULE_COUNT; i++) {
        grammar[i].head = -(mt_rand() % NUM_NON_TERMINALS + 1); // Random non-terminal
        grammar[i].body_len = (mt_rand() % MAX_RULE_BODY) + 1;

        for (int j = 0; j < grammar[i].body_len; j++) {
            // 50/50 chance of a terminal or non-terminal
            if (mt_rand() % 2 == 0) {
                grammar[i].body[j] = (mt_rand() % NUM_TERMINALS) + 1; // Terminals are > 0
            } else {
                grammar[i].body[j] = -(mt_rand() % NUM_NON_TERMINALS + 1); // Non-terminals are < 0
            }
        }
        // Ensure at least one rule can produce the start symbol to avoid trivial non-parses
        if (i == 1) {
            grammar[i].head = START_SYMBOL;
        }
    }

    // --- Generate Input Tokens ---
    input_tokens = (int*)malloc(INPUT_TOKEN_COUNT * sizeof(int));
    if (!input_tokens) { perror("malloc input_tokens"); exit(EXIT_FAILURE); }
    for (int i = 0; i < INPUT_TOKEN_COUNT; i++) {
        input_tokens[i] = (mt_rand() % NUM_TERMINALS) + 1;
    }

    // --- Initialize Chart ---
    chart = (StateSet*)malloc((INPUT_TOKEN_COUNT + 1) * sizeof(StateSet));
    if (!chart) { perror("malloc chart"); exit(EXIT_FAILURE); }

    for (int i = 0; i <= INPUT_TOKEN_COUNT; i++) {
        chart[i].states = NULL;
        chart[i].count = 0;
        chart[i].capacity = 0;
    }
}

void run_computation() {
    // Add initial state: (S' -> .S, 0) from rule 0
    State initial_state = { .rule_idx = 0, .dot_pos = 0, .origin_pos = 0 };
    add_state_to_set(initial_state, &chart[0]);

    for (int i = 0; i <= INPUT_TOKEN_COUNT; i++) {
        int current_state_idx = 0;
        while (current_state_idx < chart[i].count) {
            State s = chart[i].states[current_state_idx];
            Rule r = grammar[s.rule_idx];

            // Is the state complete?
            if (s.dot_pos >= r.body_len) {
                // --- COMPLETER ---
                // For all states in chart[s.origin_pos] that were looking for r.head
                for (int j = 0; j < chart[s.origin_pos].count; j++) {
                    State parent_state = chart[s.origin_pos].states[j];
                    Rule parent_rule = grammar[parent_state.rule_idx];
                    if (parent_state.dot_pos < parent_rule.body_len &&
                        parent_rule.body[parent_state.dot_pos] == r.head) {
                        State next_state = {
                            .rule_idx = parent_state.rule_idx,
                            .dot_pos = parent_state.dot_pos + 1,
                            .origin_pos = parent_state.origin_pos
                        };
                        add_state_to_set(next_state, &chart[i]);
                    }
                }
            } else {
                int next_symbol = r.body[s.dot_pos];
                if (next_symbol > 0) { // Terminal
                    // --- SCANNER ---
                    if (i < INPUT_TOKEN_COUNT && next_symbol == input_tokens[i]) {
                        State next_state = {
                            .rule_idx = s.rule_idx,
                            .dot_pos = s.dot_pos + 1,
                            .origin_pos = s.origin_pos
                        };
                        add_state_to_set(next_state, &chart[i+1]);
                    }
                } else { // Non-terminal
                    // --- PREDICTOR ---
                    // For all rules that start with this non-terminal
                    for (int j = 1; j <= GRAMMAR_RULE_COUNT; j++) {
                        if (grammar[j].head == next_symbol) {
                            State new_state = {
                                .rule_idx = j,
                                .dot_pos = 0,
                                .origin_pos = i
                            };
                            add_state_to_set(new_state, &chart[i]);
                        }
                    }
                }
            }
            current_state_idx++;
        }
    }

    // A simple result to prevent dead code elimination: count total states.
    long long total_states = 0;
    for (int i = 0; i <= INPUT_TOKEN_COUNT; i++) {
        total_states += chart[i].count;
    }
    final_result = total_states;
}

void cleanup() {
    free(grammar);
    free(input_tokens);
    for (int i = 0; i <= INPUT_TOKEN_COUNT; i++) {
        free(chart[i].states);
    }
    free(chart);
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
    printf("%lld\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
