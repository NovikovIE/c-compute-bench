#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <stdbool.h>

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
// --- End of MT19937 ---

// --- Benchmark Data Structures ---

typedef enum {
    TERM_TERMINAL,      // Represents a literal token (e.g., 'a', '+')
    TERM_NON_TERMINAL   // Represents a reference to another rule
} TermType;

typedef struct {
    TermType type;
    int value;          // Token ID if terminal, rule index if non-terminal
} Term;

typedef struct {
    Term* terms;
    int num_terms;
} Choice;

typedef struct {
    Choice* choices;
    int num_choices;
} Rule;

// Global state structure
struct {
    int input_token_count;
    int num_grammar_rules;
    Rule* grammar;
    int* input_tokens;
    int** memo_table; // Packrat memoization table
    long long final_result;
} state;

// --- Core Parsing Logic ---

// Recursive descent parser with memoization (Packrat parser)
// Returns number of tokens consumed, or -1 on failure.
int parse_expression(int position, int rule_index) {
    if (position >= state.input_token_count) {
        return -1; // Cannot match past the end of input
    }
    if (state.memo_table[position][rule_index] != -2) {
        return state.memo_table[position][rule_index];
    }

    Rule* rule = &state.grammar[rule_index];
    for (int i = 0; i < rule->num_choices; i++) {
        Choice* choice = &rule->choices[i];
        int current_pos = position;
        bool choice_matched = true;

        for (int j = 0; j < choice->num_terms; j++) {
            if (current_pos >= state.input_token_count) {
                choice_matched = false; // Ran out of input while matching a sequence
                break;
            }
            Term* term = &choice->terms[j];
            if (term->type == TERM_TERMINAL) {
                if (state.input_tokens[current_pos] == term->value) {
                    current_pos++;
                } else {
                    choice_matched = false;
                    break;
                }
            } else { // TERM_NON_TERMINAL
                int sub_match_len = parse_expression(current_pos, term->value);
                if (sub_match_len != -1) {
                    current_pos += sub_match_len;
                } else {
                    choice_matched = false;
                    break;
                }
            }
        }

        if (choice_matched) {
            int len_consumed = current_pos - position;
            state.memo_table[position][rule_index] = len_consumed;
            return len_consumed;
        }
    }

    state.memo_table[position][rule_index] = -1; // Failure to match this rule at this position
    return -1;
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <input_token_count> <num_grammar_rules> <seed>\n", argv[0]);
        exit(1);
    }

    state.input_token_count = atoi(argv[1]);
    state.num_grammar_rules = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);
    mt_seed(seed);

    // Generate grammar
    state.grammar = (Rule*)malloc(state.num_grammar_rules * sizeof(Rule));
    for (int i = 0; i < state.num_grammar_rules; i++) {
        state.grammar[i].num_choices = 1 + (mt_rand() % 3); // 1 to 3 choices
        state.grammar[i].choices = (Choice*)malloc(state.grammar[i].num_choices * sizeof(Choice));
        for (int j = 0; j < state.grammar[i].num_choices; j++) {
            state.grammar[i].choices[j].num_terms = 1 + (mt_rand() % 4); // 1 to 4 terms
            state.grammar[i].choices[j].terms = (Term*)malloc(state.grammar[i].choices[j].num_terms * sizeof(Term));
            for (int k = 0; k < state.grammar[i].choices[j].num_terms; k++) {
                if (mt_rand() % 2 == 0) { // 50% chance of terminal
                    state.grammar[i].choices[j].terms[k].type = TERM_TERMINAL;
                    state.grammar[i].choices[j].terms[k].value = mt_rand() % 256;
                } else {
                    state.grammar[i].choices[j].terms[k].type = TERM_NON_TERMINAL;
                    // Avoid direct left recursion (Rule -> Rule ...)
                    int target_rule;
                    do {
                        target_rule = mt_rand() % state.num_grammar_rules;
                    } while (target_rule == i);
                    state.grammar[i].choices[j].terms[k].value = target_rule;
                }
            }
        }
    }

    // Generate input token stream
    state.input_tokens = (int*)malloc(state.input_token_count * sizeof(int));
    for (int i = 0; i < state.input_token_count; i++) {
        state.input_tokens[i] = mt_rand() % 256;
    }

    // Allocate and initialize memoization table
    state.memo_table = (int**)malloc(state.input_token_count * sizeof(int*));
    for (int i = 0; i < state.input_token_count; i++) {
        state.memo_table[i] = (int*)malloc(state.num_grammar_rules * sizeof(int));
    }

    state.final_result = 0;
}

void run_computation() {
    // Initialize memoization table with a sentinel value (-2 means not computed)
    for (int i = 0; i < state.input_token_count; i++) {
        for (int j = 0; j < state.num_grammar_rules; j++) {
            state.memo_table[i][j] = -2;
        }
    }

    // Start parsing from rule 0 at position 0
    int match_length = parse_expression(0, 0);

    // Accumulate result to prevent dead code elimination
    if (match_length > 0) {
        state.final_result += match_length;
    } else {
        state.final_result += 0;
    }
}

void cleanup() {
    for (int i = 0; i < state.num_grammar_rules; i++) {
        for (int j = 0; j < state.grammar[i].num_choices; j++) {
            free(state.grammar[i].choices[j].terms);
        }
        free(state.grammar[i].choices);
    }
    free(state.grammar);

    free(state.input_tokens);

    // Free memoization table
    for (int i = 0; i < state.input_token_count; i++) {
        free(state.memo_table[i]);
    }
    free(state.memo_table);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    if (clock_gettime(CLOCK_MONOTONIC, &start) != 0) {
        perror("clock_gettime start");
        return 1;
    }

    run_computation();

    if (clock_gettime(CLOCK_MONOTONIC, &end) != 0) {
        perror("clock_gettime end");
        return 1;
    }

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%lld\n", state.final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
