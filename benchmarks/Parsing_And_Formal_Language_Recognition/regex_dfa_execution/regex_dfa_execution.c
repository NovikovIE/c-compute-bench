#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- Mersenne Twister (MT19937) Generator (Verbatim) ---
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

// --- Benchmark Globals ---
#define ALPHABET_SIZE 4

// Parameters
static long input_string_length;
static int num_dfa_states;

// Data structures
static int *input_string;
static int **transition_table;
static int *accepting_states;

// Result
static long final_result;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <input_string_length> <num_dfa_states> <seed>\n", argv[0]);
        exit(1);
    }

    input_string_length = atol(argv[1]);
    num_dfa_states = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    mt_seed(seed);

    // Allocate memory
    input_string = (int *)malloc(input_string_length * sizeof(int));
    accepting_states = (int *)malloc(num_dfa_states * sizeof(int));
    transition_table = (int **)malloc(num_dfa_states * sizeof(int *));
    if (!input_string || !accepting_states || !transition_table) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }
    for (int i = 0; i < num_dfa_states; ++i) {
        transition_table[i] = (int *)malloc(ALPHABET_SIZE * sizeof(int));
        if (!transition_table[i]) {
            fprintf(stderr, "Memory allocation failed\n");
            exit(1);
        }
    }

    // Initialize DFA: transition table and accepting states
    for (int i = 0; i < num_dfa_states; ++i) {
        // ~10% of states are accepting
        accepting_states[i] = (mt_rand() % 10 == 0) ? 1 : 0;
        for (int j = 0; j < ALPHABET_SIZE; ++j) {
            transition_table[i][j] = mt_rand() % num_dfa_states;
        }
    }

    // Generate random input string
    for (long i = 0; i < input_string_length; ++i) {
        input_string[i] = mt_rand() % ALPHABET_SIZE;
    }
}

void run_computation() {
    long accepted_count = 0;
    int current_state = 0; // Start state is always 0

    for (long i = 0; i < input_string_length; ++i) {
        int symbol = input_string[i];
        // The core of the DFA: S(t+1) = T(S(t), I(t))
        current_state = transition_table[current_state][symbol];

        // The result is an accumulation of a property observed during execution,
        // preventing dead code elimination and making the benchmark more realistic.
        // Here, we count how many times we enter an accepting state.
        if (accepting_states[current_state]) {
            accepted_count++;
        }
    }
    final_result = accepted_count;
}

void cleanup() {
    for (int i = 0; i < num_dfa_states; ++i) {
        free(transition_table[i]);
    }
    free(transition_table);
    free(accepting_states);
    free(input_string);
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
    printf("%ld\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
