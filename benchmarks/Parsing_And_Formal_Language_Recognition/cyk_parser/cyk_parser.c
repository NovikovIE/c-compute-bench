#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>

// --- Mersenne Twister (MT19937) Generator ---
// Do not modify this section
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


// --- Benchmark Data Structures ---
// Grammar rule of the form A -> B C
typedef struct {
    int A; // Index of the resulting non-terminal
    int B; // Index of the first non-terminal in the body
    int C; // Index of the second non-terminal in the body
} NonTerminalRule;

// Benchmark parameters
static int input_token_count; // n: Number of tokens in the input string
static int num_non_terminals; // r: Number of non-terminal symbols
static int grammar_rule_count_in_cnf; // g: Number of rules of the form A -> B C

// Input data
static int* input_tokens; // The input string as a sequence of token indices
static NonTerminalRule* non_terminal_rules; // The grammar rules

// DP table for CYK algorithm
// Stored flattened: P[length-1][start_pos][non_terminal]
static bool* P;

// Result
static long long final_result;

// Helper to index the flattened 3D DP table
#define P_IDX(len_idx, start_idx, nt_idx) \
    (((size_t)(len_idx) * input_token_count * num_non_terminals) + \
     ((size_t)(start_idx) * num_non_terminals) + (nt_idx))


// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <input_token_count> <num_non_terminals> <grammar_rule_count_in_cnf> <seed>\n", argv[0]);
        exit(1);
    }

    input_token_count = atoi(argv[1]);
    num_non_terminals = atoi(argv[2]);
    grammar_rule_count_in_cnf = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    if(input_token_count <= 0 || num_non_terminals <= 0 || grammar_rule_count_in_cnf <= 0) {
        fprintf(stderr, "FATAL: All parameters must be positive integers.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate memory for the input string
    input_tokens = (int*)malloc(input_token_count * sizeof(int));
    if (!input_tokens) { fprintf(stderr, "FATAL: malloc failed for input_tokens\n"); exit(1); }

    // Allocate memory for grammar rules
    non_terminal_rules = (NonTerminalRule*)malloc(grammar_rule_count_in_cnf * sizeof(NonTerminalRule));
    if (!non_terminal_rules) { fprintf(stderr, "FATAL: malloc failed for non_terminal_rules\n"); exit(1); }
    
    // Allocate memory for the CYK DP table (flattened 3D array)
    size_t p_size = (size_t)input_token_count * input_token_count * num_non_terminals;
    P = (bool*)malloc(p_size * sizeof(bool));
    if (!P) { fprintf(stderr, "FATAL: malloc failed for DP table P\n"); exit(1); }

    // Generate a random input string.
    // Each token is an integer from 0 to num_non_terminals-1.
    // This assumes the number of terminal symbols equals the number of non-terminals.
    for (int i = 0; i < input_token_count; i++) {
        input_tokens[i] = mt_rand() % num_non_terminals;
    }

    // Generate random grammar rules of the form A -> B C.
    for (int i = 0; i < grammar_rule_count_in_cnf; i++) {
        non_terminal_rules[i].A = mt_rand() % num_non_terminals;
        non_terminal_rules[i].B = mt_rand() % num_non_terminals;
        non_terminal_rules[i].C = mt_rand() % num_non_terminals;
    }
}

void run_computation() {
    // Initialize DP table to all false
    size_t p_size = (size_t)input_token_count * input_token_count * num_non_terminals;
    memset(P, 0, p_size * sizeof(bool));

    // Step 1: Handle substrings of length 1 (the base case).
    // Terminal rules are assumed to be simple: Non-terminal 'i' can produce terminal 'i'.
    for (int s = 0; s < input_token_count; s++) {
        int token_id = input_tokens[s];
        // The substring of length 1 starting at s can be derived from Non-terminal 'token_id'.
        // P[length=1][start=s][non_term=token_id] = true. `len_idx` = 1-1 = 0.
        P[P_IDX(0, s, token_id)] = true;
    }
    
    // Step 2: Fill the DP table for substrings of length 2 to n.
    // l: length of the span being considered
    for (int l = 2; l <= input_token_count; l++) {
        // s: start position of the span
        for (int s = 0; s <= input_token_count - l; s++) {
            // p: partition position of the span
            for (int p = 1; p < l; p++) {
                // For each rule A -> B C, check if the two sub-spans can be
                // generated by B and C respectively.
                for (int r = 0; r < grammar_rule_count_in_cnf; r++) {
                    int A = non_terminal_rules[r].A;
                    int B = non_terminal_rules[r].B;
                    int C = non_terminal_rules[r].C;
                    
                    // Check if P[length p][start s] contains B AND P[length l-p][start s+p] contains C.
                    bool can_form_B = P[P_IDX(p - 1, s, B)];
                    bool can_form_C = P[P_IDX(l - p - 1, s + p, C)];

                    if (can_form_B && can_form_C) {
                        // If so, the span of length l starting at s can be generated by A.
                        P[P_IDX(l - 1, s, A)] = true;
                    }
                }
            }
        }
    }

    // Step 3: Compute final result to prevent dead code elimination.
    // Sum all `true` entries in the DP table.
    final_result = 0;
    for (size_t i = 0; i < p_size; i++) {
        if (P[i]) {
            final_result++;
        }
    }
}

void cleanup() {
    free(input_tokens);
    free(non_terminal_rules);
    free(P);
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
