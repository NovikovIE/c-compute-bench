#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- Mersenne Twister (MT19937) --- (DO NOT MODIFY)
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

// --- Benchmark-specific data structures and globals ---

#define MAX_RULE_SYMBOLS 6

typedef struct {
    int* symbols;
    int length;
} Rule;

// We distinguish between non-terminals and terminals by their integer value.
// Non-terminals: [0, grammar_rule_count - 1]
// Terminals:     [grammar_rule_count, grammar_rule_count + num_terminals - 1]

static struct {
    int input_token_count;
    int grammar_rule_count;
    int num_terminals;
    int max_nesting_depth;

    int* input_tokens;
    Rule* grammar;

    size_t current_token_pos; 
    long long final_result;
} g_data;

// Forward declaration for recursive generation
void generate_stream_from_rule(int rule_id, int depth, int* token_buffer, size_t* token_idx, size_t max_tokens);

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <input_token_count> <grammar_rule_count> <max_nesting_depth> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.input_token_count = atoi(argv[1]);
    g_data.grammar_rule_count = atoi(argv[2]);
    g_data.max_nesting_depth = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    mt_seed(seed);

    // The number of unique terminal symbols will be equal to the number of rules
    g_data.num_terminals = g_data.grammar_rule_count;

    // 1. Allocate and generate the grammar
    g_data.grammar = (Rule*)malloc(g_data.grammar_rule_count * sizeof(Rule));
    if (!g_data.grammar) { perror("malloc grammar"); exit(1); }

    for (int i = 0; i < g_data.grammar_rule_count; ++i) {
        int rule_len = 1 + (mt_rand() % (MAX_RULE_SYMBOLS - 1));
        g_data.grammar[i].length = rule_len;
        g_data.grammar[i].symbols = (int*)malloc(rule_len * sizeof(int));
        if (!g_data.grammar[i].symbols) { perror("malloc rule symbols"); exit(1); }

        for (int j = 0; j < rule_len; ++j) {
             // 50/50 chance of being a terminal or non-terminal
            if (mt_rand() % 2 == 0) {
                // Non-terminal. To prevent left-recursion and guarantee termination,
                // a rule 'i' can only refer to a rule 'k' where k > i.
                int next_rule_id = i + 1 + (mt_rand() % (g_data.grammar_rule_count - i - 1));
                 if (next_rule_id >= g_data.grammar_rule_count) { // fallback to a terminal
                    g_data.grammar[i].symbols[j] = g_data.grammar_rule_count + (mt_rand() % g_data.num_terminals);
                 } else {
                    g_data.grammar[i].symbols[j] = next_rule_id;
                 }
            } else {
                // Terminal
                g_data.grammar[i].symbols[j] = g_data.grammar_rule_count + (mt_rand() % g_data.num_terminals);
            }
        }
    }

    // 2. Allocate and generate the input token stream based on the grammar
    g_data.input_tokens = (int*)malloc(g_data.input_token_count * sizeof(int));
    if (!g_data.input_tokens) { perror("malloc input_tokens"); exit(1); }

    size_t token_idx = 0;
    while (token_idx < g_data.input_token_count) {
        // Start generation from rule 0
        generate_stream_from_rule(0, 0, g_data.input_tokens, &token_idx, g_data.input_token_count);
    }
}

void generate_stream_from_rule(int rule_id, int depth, int* token_buffer, size_t* token_idx, size_t max_tokens) {
    if (depth >= g_data.max_nesting_depth || *token_idx >= max_tokens) {
        // At max depth, or if buffer is full, emit a random terminal and stop recursion.
        if (*token_idx < max_tokens) {
            token_buffer[(*token_idx)++] = g_data.grammar_rule_count + (mt_rand() % g_data.num_terminals);
        }
        return;
    }

    Rule* rule = &g_data.grammar[rule_id];
    for (int i = 0; i < rule->length && *token_idx < max_tokens; ++i) {
        int symbol = rule->symbols[i];
        if (symbol < g_data.grammar_rule_count) { // Non-terminal
            generate_stream_from_rule(symbol, depth + 1, token_buffer, token_idx, max_tokens);
        } else { // Terminal
            if (*token_idx < max_tokens) {
                 token_buffer[(*token_idx)++] = symbol;
            }
        }
    }
}

int parse_rule(int rule_id) {
    if (g_data.current_token_pos >= g_data.input_token_count) {
        return 0; // cannot match if no tokens are left
    }

    size_t saved_pos = g_data.current_token_pos;
    Rule* rule = &g_data.grammar[rule_id];

    for (int i = 0; i < rule->length; ++i) {
        int symbol = rule->symbols[i];

        if (g_data.current_token_pos >= g_data.input_token_count) {
            g_data.current_token_pos = saved_pos; // backtrack on reaching end of input mid-rule
            return 0;
        }

        if (symbol < g_data.grammar_rule_count) { // Non-terminal
            if (!parse_rule(symbol)) {
                g_data.current_token_pos = saved_pos; // backtrack on sub-rule failure
                return 0;
            }
        } else { // Terminal
            if (g_data.input_tokens[g_data.current_token_pos] == symbol) {
                g_data.current_token_pos++;
            } else {
                g_data.current_token_pos = saved_pos; // backtrack on token mismatch
                return 0;
            }
        }
    }
    return 1; // Success, entire rule matched
}

void run_computation() {
    g_data.final_result = 0;
    g_data.current_token_pos = 0;

    while (g_data.current_token_pos < g_data.input_token_count) {
        size_t initial_pos = g_data.current_token_pos;
        
        // Try to parse an expression starting from the top-level rule (rule 0).
        if (parse_rule(0)) {
            g_data.final_result++;
             // If we successfully parsed something, but consumed no tokens (e.g., empty rule), 
             // we must advance to prevent an infinite loop.
            if (g_data.current_token_pos == initial_pos) {
                 g_data.current_token_pos++;
            }
        } else {
            // On failure, employ simple error recovery: skip one token and try again.
            g_data.current_token_pos = initial_pos + 1;
        }
    }
}

void cleanup() {
    for (int i = 0; i < g_data.grammar_rule_count; ++i) {
        free(g_data.grammar[i].symbols);
    }
    free(g_data.grammar);
    free(g_data.input_tokens);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout
    printf("%lld\n", g_data.final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
