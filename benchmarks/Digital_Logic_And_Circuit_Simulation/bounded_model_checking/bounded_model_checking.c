#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- BEGIN: Mersenne Twister (DO NOT MODIFY) ---
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
// --- END: Mersenne Twister ---

#define FATAL_ERROR(msg) do { fprintf(stderr, "FATAL ERROR: %s\n", msg); exit(1); } while(0)

// Represents a single logic gate (AND, OR, XOR, NOT)
typedef enum {
    OP_AND,
    OP_OR,
    OP_XOR,
    OP_NOT,
    OP_TYPE_COUNT
} OpType;

typedef struct {
    OpType op;
    int in1;
    int in2; // Unused for NOT
} Gate;

// Global struct to hold all benchmark data
struct {
    int num_flip_flops;
    int num_gates;
    int property_complexity;
    int max_depth_k;

    // Circuit state is represented as an array of characters (booleans)
    char *state_current;
    char *state_next;

    // Definitions for the transition logic and the property to check
    Gate *transition_logic;
    Gate *property_logic;

    // Pre-allocated temporary buffers for computation
    char *transition_outputs;
    char *property_outputs;

    // Final result: number of times the property was violated
    int violation_count;
} g_data;


// Sets up the benchmark data structures and populates them with random values.
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s num_flip_flops num_gates property_complexity max_depth_k seed\n", argv[0]);
        exit(1);
    }

    g_data.num_flip_flops = atoi(argv[1]);
    g_data.num_gates = atoi(argv[2]);
    g_data.property_complexity = atoi(argv[3]);
    g_data.max_depth_k = atoi(argv[4]);
    uint32_t seed = (uint32_t)strtoul(argv[5], NULL, 10);
    
    mt_seed(seed);

    if (g_data.num_gates < g_data.num_flip_flops) {
        FATAL_ERROR("num_gates must be greater than or equal to num_flip_flops.");
    }

    // Allocate memory for all data structures
    g_data.state_current = (char*)malloc(g_data.num_flip_flops * sizeof(char));
    g_data.state_next = (char*)malloc(g_data.num_flip_flops * sizeof(char));
    g_data.transition_logic = (Gate*)malloc(g_data.num_gates * sizeof(Gate));
    g_data.property_logic = (Gate*)malloc(g_data.property_complexity * sizeof(Gate));
    g_data.transition_outputs = (char*)malloc(g_data.num_gates * sizeof(char));
    g_data.property_outputs = (char*)malloc(g_data.property_complexity * sizeof(char));

    if (!g_data.state_current || !g_data.state_next || !g_data.transition_logic || !g_data.property_logic || !g_data.transition_outputs || !g_data.property_outputs) {
         FATAL_ERROR("Memory allocation failed.");
    }

    // Initialize initial state randomly
    for (int i = 0; i < g_data.num_flip_flops; ++i) {
        g_data.state_current[i] = mt_rand() % 2;
    }

    // Generate random transition logic (a directed acyclic graph)
    for (int i = 0; i < g_data.num_gates; ++i) {
        Gate* g = &g_data.transition_logic[i];
        g->op = (OpType)(mt_rand() % OP_TYPE_COUNT);
        int num_available_inputs = g_data.num_flip_flops + i;
        g->in1 = mt_rand() % num_available_inputs;
        if (g->op != OP_NOT) {
            g->in2 = mt_rand() % num_available_inputs;
        }
    }

    // Generate random property logic
    for (int i = 0; i < g_data.property_complexity; ++i) {
        Gate* p = &g_data.property_logic[i];
        p->op = (OpType)(mt_rand() % OP_TYPE_COUNT);
        // Inputs can be any state flip-flop or any previous property gate output
        int num_available_inputs = g_data.num_flip_flops + i;
        p->in1 = mt_rand() % num_available_inputs;
        if (p->op != OP_NOT) {
            p->in2 = mt_rand() % num_available_inputs;
        }
    }
}

// Core computation: simulates the circuit for k steps.
void run_computation() {
    g_data.violation_count = 0;

    for (int k = 0; k < g_data.max_depth_k; ++k) {
        // 1. Evaluate transition logic to determine the next state
        for (int i = 0; i < g_data.num_gates; ++i) {
            Gate* g = &g_data.transition_logic[i];
            
            int in1_idx = g->in1;
            char v1 = (in1_idx < g_data.num_flip_flops) ? g_data.state_current[in1_idx] : g_data.transition_outputs[in1_idx - g_data.num_flip_flops];
            char v2 = 0;
            if (g->op != OP_NOT) {
                int in2_idx = g->in2;
                v2 = (in2_idx < g_data.num_flip_flops) ? g_data.state_current[in2_idx] : g_data.transition_outputs[in2_idx - g_data.num_flip_flops];
            }
            
            switch (g->op) {
                case OP_AND: g_data.transition_outputs[i] = v1 & v2; break;
                case OP_OR:  g_data.transition_outputs[i] = v1 | v2; break;
                case OP_XOR: g_data.transition_outputs[i] = v1 ^ v2; break;
                case OP_NOT: g_data.transition_outputs[i] = !v1;    break;
                case OP_TYPE_COUNT: break; // Should not happen
            }
        }
        // The next state is determined by the last `num_flip_flops` gate outputs
        for (int i = 0; i < g_data.num_flip_flops; ++i) {
            g_data.state_next[i] = g_data.transition_outputs[g_data.num_gates - g_data.num_flip_flops + i];
        }

        // 2. Evaluate property on the current state
        for (int i = 0; i < g_data.property_complexity; ++i) {
            Gate* p = &g_data.property_logic[i];
            
            int in1_idx = p->in1;
            char v1 = (in1_idx < g_data.num_flip_flops) ? g_data.state_current[in1_idx] : g_data.property_outputs[in1_idx - g_data.num_flip_flops];
            char v2 = 0;
            if (p->op != OP_NOT) {
                int in2_idx = p->in2;
                v2 = (in2_idx < g_data.num_flip_flops) ? g_data.state_current[in2_idx] : g_data.property_outputs[in2_idx - g_data.num_flip_flops];
            }
            
            switch (p->op) {
                 case OP_AND: g_data.property_outputs[i] = v1 & v2; break;
                 case OP_OR:  g_data.property_outputs[i] = v1 | v2; break;
                 case OP_XOR: g_data.property_outputs[i] = v1 ^ v2; break;
                 case OP_NOT: g_data.property_outputs[i] = !v1;    break;
                 case OP_TYPE_COUNT: break; // Should not happen
            }
        }
        // The property holds if the last gate's output is 1. A violation is when it's 0.
        if (g_data.property_outputs[g_data.property_complexity - 1] == 0) {
            g_data.violation_count++;
        }
        
        // 3. Advance to the next state (by swapping pointers)
        char *temp_state = g_data.state_current;
        g_data.state_current = g_data.state_next;
        g_data.state_next = temp_state;
    }
}

// Frees all memory allocated in setup_benchmark.
void cleanup() {
    free(g_data.state_current);
    free(g_data.state_next);
    free(g_data.transition_logic);
    free(g_data.property_logic);
    free(g_data.transition_outputs);
    free(g_data.property_outputs);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final result (total violations found) to stdout
    printf("%d\n", g_data.violation_count);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
