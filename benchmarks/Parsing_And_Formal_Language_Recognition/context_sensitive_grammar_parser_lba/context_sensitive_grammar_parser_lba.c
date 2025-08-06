#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

// --- START MERSENNE TWISTER (MT19937) ---
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
// --- END MERSENNE TWISTER (MT19937) ---

// Represents a simple context-sensitive rule: `left_context` `from` `right_context` -> `left_context` `to` `right_context`
typedef struct {
    char from;
    char to;
    char left_context;
    char right_context;
} Rule;

// Benchmark parameters
static int INPUT_STRING_LENGTH;
static int NUM_GRAMMAR_RULES;
static int MAX_CONFIGURATIONS;

// Global data structures
static Rule *rules;
static char **tapes_queue; // A queue of LBA tapes (configurations)

// The final result, volatile to prevent dead code elimination
static volatile long long final_result;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <input_string_length> <num_grammar_rules> <max_configurations> <seed>\n", argv[0]);
        exit(1);
    }

    INPUT_STRING_LENGTH = atoi(argv[1]);
    NUM_GRAMMAR_RULES = atoi(argv[2]);
    MAX_CONFIGURATIONS = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    mt_seed(seed);

    // Use a small alphabet to increase the probability of rule matches
    const char ALPHABET[] = "abcd";
    const int ALPHABET_SIZE = sizeof(ALPHABET) - 1;

    // Allocate and generate grammar rules
    rules = (Rule *)malloc(NUM_GRAMMAR_RULES * sizeof(Rule));
    for (int i = 0; i < NUM_GRAMMAR_RULES; i++) {
        rules[i].from = ALPHABET[mt_rand() % ALPHABET_SIZE];
        rules[i].to = ALPHABET[mt_rand() % ALPHABET_SIZE];
        rules[i].left_context = ALPHABET[mt_rand() % ALPHABET_SIZE];
        rules[i].right_context = ALPHABET[mt_rand() % ALPHABET_SIZE];
    }

    // Allocate the queue of configurations (tapes)
    tapes_queue = (char **)malloc(MAX_CONFIGURATIONS * sizeof(char *));
    for (int i = 0; i < MAX_CONFIGURATIONS; i++) {
        // +1 for null terminator
        tapes_queue[i] = (char *)malloc((INPUT_STRING_LENGTH + 1) * sizeof(char));
    }

    // Create the initial configuration (tape)
    for (int i = 0; i < INPUT_STRING_LENGTH; i++) {
        tapes_queue[0][i] = ALPHABET[mt_rand() % ALPHABET_SIZE];
    }
    tapes_queue[0][INPUT_STRING_LENGTH] = '\0';
}

void run_computation() {
    int head = 0; // The configuration we are currently processing
    int tail = 1; // The next available slot for a new configuration
    long long result_accumulator = 0;

    // The main simulation loop of the Linear Bounded Automaton
    // Process configurations from the queue until it's empty or we run out of space
    while (head < tail && tail < MAX_CONFIGURATIONS) {
        char *current_tape = tapes_queue[head];

        // Try to apply each rule at each position on the tape
        for (int i = 0; i < INPUT_STRING_LENGTH; i++) {
            for (int r = 0; r < NUM_GRAMMAR_RULES; r++) {
                if (current_tape[i] != rules[r].from) {
                    continue;
                }

                // Check for context
                char current_left = (i > 0) ? current_tape[i - 1] : -1; // Use -1 for no context
                char current_right = (i < INPUT_STRING_LENGTH - 1) ? current_tape[i + 1] : -1;

                if (current_left == rules[r].left_context && current_right == rules[r].right_context) {
                    // Rule matches. Generate a new configuration.
                    if (tail < MAX_CONFIGURATIONS) {
                        char *new_tape = tapes_queue[tail];
                        memcpy(new_tape, current_tape, INPUT_STRING_LENGTH + 1);
                        new_tape[i] = rules[r].to;

                        // Accumulate a value to prevent optimization
                        result_accumulator += new_tape[i];
                        
                        tail++;
                    } else {
                        // Configuration space is full, stop generating new ones
                        goto end_simulation;
                    }
                }
            }
        }
        head++;
    }

end_simulation:
    final_result = result_accumulator;
}

void cleanup() {
    free(rules);
    for (int i = 0; i < MAX_CONFIGURATIONS; i++) {
        free(tapes_queue[i]);
    }
    free(tapes_queue);
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
