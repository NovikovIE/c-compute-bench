/**
 * @file regex_nfa_simulation.c
 * @brief Benchmark for simulating a Non-deterministic Finite Automaton (NFA).
 * 
 * @details This program simulates the process of recognizing a string using an NFA,
 * a core concept in formal language theory and compiler design (lexical analysis).
 * The NFA's transition table and the input string are generated randomly.
 * The simulation proceeds by tracking the set of all possible current states
 * as it consumes the input string one character at a time.
 *
 * The complexity of the simulation is O(input_string_length * num_nfa_states^2)
 * for this dense-matrix implementation.
 *
 * @param input_string_length The length of the random input string to be processed.
 * @param num_nfa_states The total number of states in the NFA.
 * @param seed The seed for the random number generator.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <string.h>

// --- START: Mersenne Twister (MT19937) --- Do Not Modify ---
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
// --- END: Mersenne Twister (MT19937) ---

#define ALPHABET_SIZE 2 // Using a simple alphabet {'a', 'b'}
#define TRANSITION_DENSITY_PPT 2 // Parts-per-thousand for transition probability (0.2%)
#define ACCEPT_STATE_PROB 10 // Probability (in %) for a state to be an accepting state

// Global struct to hold all benchmark data
typedef struct {
    int input_string_length;
    int num_nfa_states;

    char *input_string;

    // NFA Representation
    // transitions[from_state][char_index][to_state] is a boolean (char)
    // char_index 0 for 'a', 1 for 'b'
    char*** transitions;
    char* accept_states; // boolean (char) array
    
    int final_result; // To store the result from run_computation
} BenchmarkData;

static BenchmarkData* g_data;

/**
 * @brief Parses command line arguments, allocates memory, and generates NFA and input string.
 */
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <input_string_length> <num_nfa_states> <seed>\n", argv[0]);
        exit(1);
    }

    g_data = (BenchmarkData*)malloc(sizeof(BenchmarkData));
    if (!g_data) {
        perror("Failed to allocate memory for benchmark data");
        exit(1);
    }

    g_data->input_string_length = atoi(argv[1]);
    g_data->num_nfa_states = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    mt_seed(seed);

    // Allocate and generate the input string
    g_data->input_string = (char*)malloc((g_data->input_string_length + 1) * sizeof(char));
    for (int i = 0; i < g_data->input_string_length; ++i) {
        g_data->input_string[i] = 'a' + (mt_rand() % ALPHABET_SIZE);
    }
    g_data->input_string[g_data->input_string_length] = '\0';

    // Allocate and generate the NFA's accepting states
    g_data->accept_states = (char*)malloc(g_data->num_nfa_states * sizeof(char));
    for (int i = 0; i < g_data->num_nfa_states; ++i) {
        g_data->accept_states[i] = (mt_rand() % 100 < ACCEPT_STATE_PROB);
    }

    // Allocate and generate the NFA transition table (dense matrix representation)
    g_data->transitions = (char***)malloc(g_data->num_nfa_states * sizeof(char**));
    for (int i = 0; i < g_data->num_nfa_states; ++i) {
        g_data->transitions[i] = (char**)malloc(ALPHABET_SIZE * sizeof(char*));
        for (int j = 0; j < ALPHABET_SIZE; ++j) {
            g_data->transitions[i][j] = (char*)malloc(g_data->num_nfa_states * sizeof(char));
            for (int k = 0; k < g_data->num_nfa_states; ++k) {
                g_data->transitions[i][j][k] = (mt_rand() % 1000 < TRANSITION_DENSITY_PPT);
            }
        }
    }
}

/**
 * @brief Executes the NFA simulation on the generated data.
 */
void run_computation() {
    char* current_states = (char*)calloc(g_data->num_nfa_states, sizeof(char));
    char* next_states = (char*)calloc(g_data->num_nfa_states, sizeof(char));
    if (!current_states || !next_states) {
        perror("Failed to allocate memory for state arrays");
        exit(1);
    }

    // Start state is always state 0
    current_states[0] = 1;

    // Process each character in the input string
    for (int i = 0; i < g_data->input_string_length; ++i) {
        memset(next_states, 0, g_data->num_nfa_states * sizeof(char));
        int char_idx = g_data->input_string[i] - 'a';

        // For each currently possible state, find all possible next states
        for (int s = 0; s < g_data->num_nfa_states; ++s) {
            if (current_states[s]) {
                // Check for transitions to all other states
                for (int s_next = 0; s_next < g_data->num_nfa_states; ++s_next) {
                    if (g_data->transitions[s][char_idx][s_next]) {
                        next_states[s_next] = 1;
                    }
                }
            }
        }

        // The set of next states becomes the new set of current states
        memcpy(current_states, next_states, g_data->num_nfa_states * sizeof(char));
    }

    // After processing the whole string, count how many final states are accept states
    int accepted_count = 0;
    for (int s = 0; s < g_data->num_nfa_states; ++s) {
        if (current_states[s] && g_data->accept_states[s]) {
            accepted_count++;
        }
    }
    g_data->final_result = accepted_count;

    free(current_states);
    free(next_states);
}

/**
 * @brief Frees all memory allocated in setup_benchmark.
 */
void cleanup() {
    for (int i = 0; i < g_data->num_nfa_states; ++i) {
        for (int j = 0; j < ALPHABET_SIZE; ++j) {
            free(g_data->transitions[i][j]);
        }
        free(g_data->transitions[i]);
    }
    free(g_data->transitions);
    free(g_data->accept_states);
    free(g_data->input_string);
    free(g_data);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout to prevent dead code elimination
    printf("%d\n", g_data->final_result);

    cleanup();

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
