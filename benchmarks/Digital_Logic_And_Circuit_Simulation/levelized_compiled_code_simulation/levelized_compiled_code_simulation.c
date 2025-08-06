#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

//
// --- Mersenne Twister (MT19937) Generator (DO NOT MODIFY) ---
//
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

//
// --- Benchmark: levelized_compiled_code_simulation ---
//

// Gate operation types
typedef enum {
    AND,
    OR,
    NOT,
    XOR
} GateOp;

// A single logic gate in the circuit
typedef struct {
    GateOp op;
    // Indices of inputs. A negative index -(i+1) refers to primary input i.
    // A non-negative index j refers to the output of gate j.
    int input1_idx;
    int input2_idx; // Unused for NOT gates
} Gate;

// --- Global Data Structures ---

// Benchmark parameters
static int num_gates;
static int num_input_vectors;
static int num_primary_inputs;

// Data arrays
static Gate* circuit = NULL;
static uint8_t** input_vectors = NULL;
static uint8_t* gate_values = NULL;

// Result accumulator to prevent dead-code elimination
static unsigned long long final_result = 0;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_gates> <num_input_vectors> <num_primary_inputs> <seed>\n", argv[0]);
        exit(1);
    }

    num_gates = atoi(argv[1]);
    num_input_vectors = atoi(argv[2]);
    num_primary_inputs = atoi(argv[3]);
    uint32_t seed = atoi(argv[4]);
    
    mt_seed(seed);

    // Allocate memory
    circuit = (Gate*)malloc(num_gates * sizeof(Gate));
    gate_values = (uint8_t*)malloc(num_gates * sizeof(uint8_t));

    // Allocate input vectors as a contiguous block of memory for better cache performance
    // and then set up row pointers.
    uint8_t* vectors_data = (uint8_t*)malloc((size_t)num_input_vectors * num_primary_inputs * sizeof(uint8_t));
    input_vectors = (uint8_t**)malloc(num_input_vectors * sizeof(uint8_t*));
    for (int i = 0; i < num_input_vectors; ++i) {
        input_vectors[i] = &vectors_data[i * num_primary_inputs];
    }

    if (!circuit || !gate_values || !vectors_data || !input_vectors) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Generate a random levelized circuit (DAG)
    for (int i = 0; i < num_gates; ++i) {
        circuit[i].op = (GateOp)(mt_rand() % 4);
        int num_available_sources = num_primary_inputs + i;

        // Choose input 1. Can be a primary input or output of a previous gate.
        int source1 = mt_rand() % num_available_sources;
        circuit[i].input1_idx = (source1 < num_primary_inputs) ? -(source1 + 1) : (source1 - num_primary_inputs);
        
        if (circuit[i].op != NOT) {
            // Choose input 2 for binary gates
            int source2 = mt_rand() % num_available_sources;
            circuit[i].input2_idx = (source2 < num_primary_inputs) ? -(source2 + 1) : (source2 - num_primary_inputs);
        } else {
            circuit[i].input2_idx = 0; // Unused
        }
    }

    // Generate random input vectors
    for (int i = 0; i < num_input_vectors; ++i) {
        for (int j = 0; j < num_primary_inputs; ++j) {
            input_vectors[i][j] = mt_rand() % 2;
        }
    }
}

void run_computation() {
    final_result = 0;
    for (int i = 0; i < num_input_vectors; ++i) {
        uint8_t* current_vector = input_vectors[i];

        // Simulate each gate in levelized order
        for (int g = 0; g < num_gates; ++g) {
            // Fetch input values
            int idx1 = circuit[g].input1_idx;
            uint8_t val1 = (idx1 < 0) ? current_vector[-idx1 - 1] : gate_values[idx1];

            uint8_t result;
            switch (circuit[g].op) {
                case AND: {
                    int idx2 = circuit[g].input2_idx;
                    uint8_t val2 = (idx2 < 0) ? current_vector[-idx2 - 1] : gate_values[idx2];
                    result = val1 & val2;
                    break;
                }
                case OR: {
                    int idx2 = circuit[g].input2_idx;
                    uint8_t val2 = (idx2 < 0) ? current_vector[-idx2 - 1] : gate_values[idx2];
                    result = val1 | val2;
                    break;
                }
                case XOR: {
                    int idx2 = circuit[g].input2_idx;
                    uint8_t val2 = (idx2 < 0) ? current_vector[-idx2 - 1] : gate_values[idx2];
                    result = val1 ^ val2;
                    break;
                }
                case NOT: {
                    result = !val1;
                    break;
                }
            }
            gate_values[g] = result;
        }

        // Accumulate the output of the last gate to prevent D.C.E
        if (num_gates > 0) {
            final_result += gate_values[num_gates - 1];
        }
    }
}

void cleanup() {
    if (input_vectors) {
        free(input_vectors[0]); // Free the contiguous data block
        free(input_vectors);    // Free the row pointers
    }
    free(circuit);
    free(gate_values);
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
    printf("%llu\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
