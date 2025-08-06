#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- BEGIN MERSENNE TWISTER (DO NOT MODIFY) ---
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
// --- END MERSENNE TWISTER ---

// --- BENCHMARK SPECIFIC CODE ---

// Represents the type of a logic gate.
typedef enum {
    INPUT,
    NOT,
    AND,
    OR,
    XOR
} GateType;

// Represents a single logic gate.
typedef struct {
    GateType type;
    int input1;
    int input2;       // Not used by INPUT or NOT gates
    int output_value; // Current value (0 or 1)
} Gate;

// Represents a complete combinational circuit.
typedef struct {
    Gate* gates;
    int num_gates;
    int num_inputs;
} Circuit;

// Benchmark parameters
int P_NUM_INPUTS;
int P_NUM_GATES_A;
int P_NUM_GATES_B;

// Globally accessible benchmark data
Circuit* circuit_a = NULL;
Circuit* circuit_b = NULL;
int mismatch_count = 0; // The final result

// Helper to generate a random circuit.
void generate_circuit(Circuit* c, int num_gates, int num_inputs) {
    c->gates = (Gate*)malloc(sizeof(Gate) * num_gates);
    if (!c->gates) {
        fprintf(stderr, "Failed to allocate memory for gates.\n");
        exit(1);
    }
    c->num_gates = num_gates;
    c->num_inputs = num_inputs;

    // The first `num_inputs` gates are primary inputs.
    for (int i = 0; i < num_inputs; ++i) {
        c->gates[i].type = INPUT;
    }

    // The remaining gates are randomly generated logic gates.
    // We ensure inputs to gate `i` come from gates with index < `i`,
    // creating a valid Directed Acyclic Graph (DAG).
    for (int i = num_inputs; i < num_gates; ++i) {
        // Type can be NOT, AND, OR, XOR.
        c->gates[i].type = (GateType)(1 + (mt_rand() % 4));
        c->gates[i].input1 = mt_rand() % i;
        if (c->gates[i].type != NOT) { // Binary gates need a second input
            c->gates[i].input2 = mt_rand() % i;
        }
    }
}

// Helper to evaluate a circuit's output for the current input values.
int evaluate_circuit(Circuit* c) {
    for (int i = c->num_inputs; i < c->num_gates; ++i) {
        Gate* g = &c->gates[i];
        int val1 = c->gates[g->input1].output_value;
        
        switch(g->type) {
            case NOT:
                g->output_value = !val1;
                break;
            case AND:
                g->output_value = val1 && c->gates[g->input2].output_value;
                break;
            case OR:
                g->output_value = val1 || c->gates[g->input2].output_value;
                break;
            case XOR:
                g->output_value = val1 ^ c->gates[g->input2].output_value;
                break;
            case INPUT:
                // INPUT gates are set externally, do nothing here.
                break;
        }
    }
    // The circuit's output is the output of the very last gate.
    return c->gates[c->num_gates - 1].output_value;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_inputs> <num_gates_a> <num_gates_b> <seed>\n", argv[0]);
        exit(1);
    }

    P_NUM_INPUTS = atoi(argv[1]);
    P_NUM_GATES_A = atoi(argv[2]);
    P_NUM_GATES_B = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    mt_seed(seed);

    if (P_NUM_GATES_A < P_NUM_INPUTS || P_NUM_GATES_B < P_NUM_INPUTS) {
        fprintf(stderr, "Error: Number of gates must be >= number of inputs.\n");
        exit(1);
    }
    if (P_NUM_INPUTS > 30) { // Practically limited to avoid extreme runtimes
        fprintf(stderr, "Error: num_inputs > 30 is not recommended.\n");
        exit(1);
    }
    
    circuit_a = (Circuit*)malloc(sizeof(Circuit));
    circuit_b = (Circuit*)malloc(sizeof(Circuit));
    if (!circuit_a || !circuit_b) {
        fprintf(stderr, "Failed to allocate memory for circuits.\n");
        exit(1);
    }

    generate_circuit(circuit_a, P_NUM_GATES_A, P_NUM_INPUTS);
    generate_circuit(circuit_b, P_NUM_GATES_B, P_NUM_INPUTS);
}

void run_computation() {
    uint64_t num_combinations = 1ULL << P_NUM_INPUTS;
    mismatch_count = 0;

    // Iterate through all 2^N possible input combinations.
    for (uint64_t i = 0; i < num_combinations; ++i) {
        // Set the output values of the input gates for this combination.
        for (int j = 0; j < P_NUM_INPUTS; ++j) {
            int input_val = (i >> j) & 1;
            circuit_a->gates[j].output_value = input_val;
            circuit_b->gates[j].output_value = input_val;
        }

        // Evaluate both circuits with the current inputs.
        int output_a = evaluate_circuit(circuit_a);
        int output_b = evaluate_circuit(circuit_b);

        // If the outputs differ, they are not equivalent for this input.
        if (output_a != output_b) {
            mismatch_count++;
        }
    }
}

void cleanup() {
    if (circuit_a) {
        free(circuit_a->gates);
        free(circuit_a);
    }
    if (circuit_b) {
        free(circuit_b->gates);
        free(circuit_b);
    }
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final result (total number of mismatches) to stdout.
    printf("%d\n", mismatch_count);

    // Print the time taken to stderr.
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
