#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// Mersenne Twister (MT19937) Generator
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
// End of Mersenne Twister

// Benchmark parameters
static int NUM_STATE_VARS;
static int TRANSITION_LOGIC_COMPLEXITY;
static int UNROLLING_BOUND_K;

// Data structures
// Represents a simple 3-input logic gate that determines a part of the next state
typedef struct {
    int input1;
    int input2;
    int input3;
} LogicGate;

// Represents the transition system logic T(s, s') as a collection of gates
static LogicGate **transition_logic;

// Represents the unrolled states s_0, s_1, ..., s_k
// Using char as a boolean (0 or 1) for memory efficiency
static char **states;

// Final result to prevent dead code elimination through aggregation
static unsigned int final_result;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_state_variables> <transition_logic_complexity> <unrolling_bound_k> <seed>\n", argv[0]);
        exit(1);
    }

    NUM_STATE_VARS = atoi(argv[1]);
    TRANSITION_LOGIC_COMPLEXITY = atoi(argv[2]);
    UNROLLING_BOUND_K = atoi(argv[3]);
    uint32_t seed = atoi(argv[4]);

    mt_seed(seed);

    // Allocate memory for the unrolled states
    states = (char **)malloc((UNROLLING_BOUND_K + 1) * sizeof(char *));
    for (int k = 0; k <= UNROLLING_BOUND_K; ++k) {
        states[k] = (char *)malloc(NUM_STATE_VARS * sizeof(char));
    }

    // Allocate memory for the transition logic
    transition_logic = (LogicGate **)malloc(NUM_STATE_VARS * sizeof(LogicGate *));
    for (int v = 0; v < NUM_STATE_VARS; ++v) {
        transition_logic[v] = (LogicGate *)malloc(TRANSITION_LOGIC_COMPLEXITY * sizeof(LogicGate));
    }

    // Initialize the initial state (s_0) with random boolean values
    for (int v = 0; v < NUM_STATE_VARS; ++v) {
        states[0][v] = mt_rand() % 2;
    }

    // Initialize the transition logic with random connections to previous state variables
    for (int v = 0; v < NUM_STATE_VARS; ++v) {
        for (int c = 0; c < TRANSITION_LOGIC_COMPLEXITY; ++c) {
            transition_logic[v][c].input1 = mt_rand() % NUM_STATE_VARS;
            transition_logic[v][c].input2 = mt_rand() % NUM_STATE_VARS;
            transition_logic[v][c].input3 = mt_rand() % NUM_STATE_VARS;
        }
    }
}

void run_computation() {
    // This loop simulates the core task of a Bounded Model Checker: unrolling the system's
    // transition relation for 'k' steps to create a large propositional formula.
    for (int k = 0; k < UNROLLING_BOUND_K; ++k) {
        for (int v = 0; v < NUM_STATE_VARS; ++v) {
            char next_state_val = 0;
            // The value of state variable 'v' at time 'k+1' is determined by
            // applying its associated logic gates to the state at time 'k'.
            for (int c = 0; c < TRANSITION_LOGIC_COMPLEXITY; ++c) {
                LogicGate gate = transition_logic[v][c];
                char in1 = states[k][gate.input1];
                char in2 = states[k][gate.input2];
                char in3 = states[k][gate.input3];
                
                // An arbitrary boolean function to simulate complex transition logic
                char gate_output = (in1 & !in2) | (in2 ^ in3);
                
                // Combine outputs of all gates for this variable (e.g., with XOR)
                next_state_val ^= gate_output;
            }
            states[k + 1][v] = next_state_val;
        }
    }

    // A final pass to aggregate results, preventing the compiler from optimizing away the state computations.
    // This simulates a final property check that inspects the entire execution trace.
    unsigned int checksum = 0;
    const unsigned int prime_multiplier = 31;
    for (int k = 1; k <= UNROLLING_BOUND_K; ++k) {
        for (int v = 0; v < NUM_STATE_VARS; ++v) {
            checksum = checksum * prime_multiplier + states[k][v];
        }
    }
    final_result = checksum;
}

void cleanup() {
    for (int k = 0; k <= UNROLLING_BOUND_K; ++k) {
        free(states[k]);
    }
    free(states);

    for (int v = 0; v < NUM_STATE_VARS; ++v) {
        free(transition_logic[v]);
    }
    free(transition_logic);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final aggregated result to stdout
    printf("%u\n", final_result);

    // Print timing information to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
