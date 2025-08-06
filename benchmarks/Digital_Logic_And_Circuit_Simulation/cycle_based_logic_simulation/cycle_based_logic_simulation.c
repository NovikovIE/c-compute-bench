#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <string.h>


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

// --- Benchmark Globals ---

// Parameters
int num_gates;
int num_flip_flops;
int num_clock_cycles;

// Derived
int total_nodes;

// Logic gate definition
typedef enum { AND, OR, XOR } GateType;
typedef struct {
    GateType type;
    int input1_idx;
    int input2_idx;
} Gate;

// Circuit data structures
Gate* gates;
int* ff_d_inputs;      // Input index for each flip-flop
char* current_state;   // Boolean state (0 or 1)
char* next_state;      // Next state calculated in a cycle

// Final result accumulator
long long final_result = 0;


void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_gates> <num_flip_flops> <num_clock_cycles> <seed>\n", argv[0]);
        exit(1);
    }

    num_gates = atoi(argv[1]);
    num_flip_flops = atoi(argv[2]);
    num_clock_cycles = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    
    mt_seed(seed);

    total_nodes = num_gates + num_flip_flops;

    gates = (Gate*)malloc(num_gates * sizeof(Gate));
    ff_d_inputs = (int*)malloc(num_flip_flops * sizeof(int));
    current_state = (char*)malloc(total_nodes * sizeof(char));
    next_state = (char*)malloc(total_nodes * sizeof(char));

    if (!gates || !ff_d_inputs || !current_state || !next_state) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(1);
    }

    // Initialize circuit topology randomly
    // Gate outputs are nodes 0 to num_gates-1
    for (int i = 0; i < num_gates; ++i) {
        gates[i].type = (GateType)(mt_rand() % 3); // AND, OR, or XOR
        gates[i].input1_idx = mt_rand() % total_nodes;
        gates[i].input2_idx = mt_rand() % total_nodes;
    }

    // Flip-flop outputs are nodes num_gates to total_nodes-1
    for (int i = 0; i < num_flip_flops; ++i) {
        ff_d_inputs[i] = mt_rand() % total_nodes;
    }

    // Initialize the circuit's starting state randomly
    for (int i = 0; i < total_nodes; ++i) {
        current_state[i] = mt_rand() % 2;
    }
}

void run_computation() {
    for (int cycle = 0; cycle < num_clock_cycles; ++cycle) {
        // 1. Evaluate combinational logic (gates)
        for (int i = 0; i < num_gates; ++i) {
            char in1 = current_state[gates[i].input1_idx];
            char in2 = current_state[gates[i].input2_idx];
            char out;
            switch (gates[i].type) {
                case AND: out = in1 & in2; break;
                case OR:  out = in1 | in2; break;
                case XOR: out = in1 ^ in2; break;
            }
            // The output of gate `i` defines the next state of node `i`.
            next_state[i] = out;
        }

        // 2. Evaluate inputs to sequential logic (flip-flops)
        for (int i = 0; i < num_flip_flops; ++i) {
            char d_in = current_state[ff_d_inputs[i]];
            // The D input determines the next state of the flip-flop's output node.
            // Node index for FF `i` is `num_gates + i`.
            next_state[num_gates + i] = d_in;
        }

        // 3. "Clock tick": update the state of all nodes for the next cycle.
        memcpy(current_state, next_state, total_nodes * sizeof(char));
    }

    // Calculate a final result to prevent dead code elimination.
    // The sum of all node states in the final cycle.
    long long sum = 0;
    for (int i = 0; i < total_nodes; ++i) {
        sum += current_state[i];
    }
    final_result = sum;
}

void cleanup() {
    free(gates);
    free(ff_d_inputs);
    free(current_state);
    free(next_state);
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
