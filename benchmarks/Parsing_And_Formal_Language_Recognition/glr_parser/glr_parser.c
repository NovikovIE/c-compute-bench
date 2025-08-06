#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

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

// Defines for the benchmark logic
#define ACTION_SHIFT 0
#define ACTION_REDUCE 1
#define ACTION_NONE 2
#define MAX_RULE_BODY_LEN 8

// Data structures for the GLR parser simulation
typedef struct {
    int type;  // ACTION_SHIFT, ACTION_REDUCE, or ACTION_NONE
    int value; // target state for SHIFT, rule index for REDUCE
} Action;

typedef struct {
    int head_symbol;
    int body_length;
} Rule;

typedef struct {
    int state;
} ParserProcess;

// Global structure to hold all benchmark data
typedef struct {
    int input_token_count;
    int grammar_rule_count;
    int num_states;
    int num_symbols;
    int num_terminals;
    int num_non_terminals;
    int max_processes;

    int* input_tokens;
    Rule* rules;
    Action** action_table;
    int** goto_table;

    ParserProcess* active_processes;
    ParserProcess* next_processes;
    int active_process_count;

    long long final_result;
} BenchmarkData;

static BenchmarkData* data;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <input_token_count> <grammar_rule_count> <num_states> <num_symbols> <seed>\n", argv[0]);
        exit(1);
    }

    data = (BenchmarkData*)malloc(sizeof(BenchmarkData));
    if (!data) {
        perror("Failed to allocate memory for BenchmarkData");
        exit(1);
    }

    data->input_token_count = atoi(argv[1]);
    data->grammar_rule_count = atoi(argv[2]);
    data->num_states = atoi(argv[3]);
    data->num_symbols = atoi(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    mt_seed(seed);

    data->num_terminals = data->num_symbols / 2;
    data->num_non_terminals = data->num_symbols - data->num_terminals;
    data->max_processes = 2 * data->num_states; // A heuristic cap

    // Allocate memory for all data structures
    data->input_tokens = (int*)malloc(data->input_token_count * sizeof(int));
    data->rules = (Rule*)malloc(data->grammar_rule_count * sizeof(Rule));
    data->active_processes = (ParserProcess*)malloc(data->max_processes * sizeof(ParserProcess));
    data->next_processes = (ParserProcess*)malloc(data->max_processes * sizeof(ParserProcess));

    data->action_table = (Action**)malloc(data->num_states * sizeof(Action*));
    for (int i = 0; i < data->num_states; ++i) {
        data->action_table[i] = (Action*)malloc(data->num_terminals * sizeof(Action));
    }

    data->goto_table = (int**)malloc(data->num_states * sizeof(int*));
        for (int i = 0; i < data->num_states; ++i) {
        data->goto_table[i] = (int*)malloc(data->num_non_terminals * sizeof(int));
    }

    // Generate random input tokens (terminals only)
    for (int i = 0; i < data->input_token_count; ++i) {
        data->input_tokens[i] = mt_rand() % data->num_terminals;
    }

    // Generate random grammar rules
    for (int i = 0; i < data->grammar_rule_count; ++i) {
        data->rules[i].head_symbol = (mt_rand() % data->num_non_terminals) + data->num_terminals;
        data->rules[i].body_length = (mt_rand() % MAX_RULE_BODY_LEN) + 1;
    }

    // Generate random ACTION table
    for (int i = 0; i < data->num_states; ++i) {
        for (int j = 0; j < data->num_terminals; ++j) {
            int p = mt_rand() % 100;
            if (p < 48) { // SHIFT
                data->action_table[i][j].type = ACTION_SHIFT;
                data->action_table[i][j].value = mt_rand() % data->num_states;
            } else if (p < 96) { // REDUCE
                data->action_table[i][j].type = ACTION_REDUCE;
                data->action_table[i][j].value = mt_rand() % data->grammar_rule_count;
            } else { // NOP
                data->action_table[i][j].type = ACTION_NONE;
            }
        }
    }

    // Generate random GOTO table
    for (int i = 0; i < data->num_states; ++i) {
        for (int j = 0; j < data->num_non_terminals; ++j) {
            data->goto_table[i][j] = mt_rand() % data->num_states;
        }
    }

    data->final_result = 0;
}

void run_computation() {
    // Start with one process in state 0
    data->active_processes[0].state = 0;
    data->active_process_count = 1;

    for (int i = 0; i < data->input_token_count; ++i) {
        int token = data->input_tokens[i];
        int next_process_count = 0;

        for (int j = 0; j < data->active_process_count; ++j) {
            int current_state = data->active_processes[j].state;
            Action action = data->action_table[current_state][token];

            switch (action.type) {
                case ACTION_SHIFT:
                    if (next_process_count < data->max_processes) {
                        data->next_processes[next_process_count++].state = action.value;
                    }
                    break;
                
                case ACTION_REDUCE:
                    {
                        Rule rule = data->rules[action.value];
                        // Simulate reduction work (traverse back) and GOTO
                        // In a real GLR parser, this involves complex GSS traversal.
                        // Here, we simulate the computational cost.
                        int state_after_pop = current_state;
                        for (int k = 0; k < rule.body_length; k++) {
                            state_after_pop = data->goto_table[state_after_pop][k % data->num_non_terminals];
                        }

                        int non_terminal_idx = rule.head_symbol - data->num_terminals;
                        int new_state = data->goto_table[state_after_pop][non_terminal_idx];

                        if (next_process_count < data->max_processes) {
                            data->next_processes[next_process_count++].state = new_state;
                        }
                    }
                    break;

                case ACTION_NONE:
                    // This process dies
                    break;
            }
        }

        // Swap active and next process lists
        ParserProcess* temp_ptr = data->active_processes;
        data->active_processes = data->next_processes;
        data->next_processes = temp_ptr;
        data->active_process_count = next_process_count;

        if (data->active_process_count == 0) {
            // Parsing failed, but for benchmark purposes, re-initialize to avoid early exit
            data->active_processes[0].state = 0;
            data->active_process_count = 1;
        }
    }

    // Accumulate a final result to prevent dead code elimination
    long long sum = 0;
    for (int i = 0; i < data->active_process_count; ++i) {
        sum += data->active_processes[i].state;
    }
    data->final_result = sum;
}

void cleanup() {
    for (int i = 0; i < data->num_states; ++i) {
        free(data->action_table[i]);
        free(data->goto_table[i]);
    }
    free(data->action_table);
    free(data->goto_table);

    free(data->input_tokens);
    free(data->rules);
    free(data->active_processes);
    free(data->next_processes);
    free(data);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    long long final_result = data->final_result;
    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%lld\n", final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
