#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

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
// --- END MERSENNE TWISTER ---

// --- BENCHMARK DATA STRUCTURES ---

// A 'literal' represents a wire in the circuit. It comprises a variable index and an inversion flag.
// bit 0:      is_inverted (1 = true, 0 = false)
// bits 1-31:  variable_index (identifies a primary input or a gate output)
typedef uint32_t Literal;

// A 2-input AND gate.
typedef struct {
    Literal input1;
    Literal input2;
} Gate;

// An And-Inverter Graph (AIG) representing a logic circuit.
// The graph is topologically sorted by construction.
typedef struct {
    Gate* gates;
    uint32_t num_pis;       // Number of primary inputs
    uint32_t num_gates;     // Number of AND gates
    Literal output_literal; // The single output of the circuit
} Circuit;

// --- GLOBAL DATA ---

// Benchmark parameters
static uint32_t NUM_PRIMARY_INPUTS;
static uint32_t NUM_GATES_CIRCUIT_A;
static uint32_t NUM_GATES_CIRCUIT_B;
static uint32_t NUM_SIMULATION_VECTORS;

// Global data structures for the circuits and simulation vectors
static Circuit* g_circuit_a = NULL;
static Circuit* g_circuit_b = NULL;
static uint64_t* g_simulation_vectors = NULL;

// --- HELPER FUNCTIONS ---

// Helper to create a single random circuit (AIG)
Circuit* create_random_circuit(uint32_t num_pis, uint32_t num_gates) {
    Circuit* c = (Circuit*)malloc(sizeof(Circuit));
    if (!c) return NULL;

    c->num_pis = num_pis;
    c->num_gates = num_gates;
    c->gates = (Gate*)malloc(num_gates * sizeof(Gate));
    if (!c->gates) {
        free(c);
        return NULL;
    }

    // Generate gates in topological order
    for (uint32_t i = 0; i < num_gates; ++i) {
        uint32_t max_var_index = num_pis + i;
        
        // Pick two inputs from PIs or previous gates
        uint32_t idx1 = mt_rand() % max_var_index;
        uint32_t idx2 = mt_rand() % max_var_index;
        
        // Randomly invert inputs
        uint8_t inv1 = mt_rand() & 1;
        uint8_t inv2 = mt_rand() & 1;

        c->gates[i].input1 = (idx1 << 1) | inv1;
        c->gates[i].input2 = (idx2 << 1) | inv2;
    }

    // Set the final output to be the output of the last gate, possibly inverted
    uint32_t out_idx = num_pis + num_gates - 1;
    uint8_t out_inv = mt_rand() & 1;
    c->output_literal = (out_idx << 1) | out_inv;

    return c;
}

// Helper to evaluate a circuit for a given input vector
uint8_t evaluate_circuit(const Circuit* c, uint64_t input_vector, uint8_t* gate_values) {
    for (uint32_t i = 0; i < c->num_gates; ++i) {
        Gate g = c->gates[i];

        // Evaluate first input literal
        uint32_t var_idx1 = g.input1 >> 1;
        uint8_t is_inv1 = g.input1 & 1;
        uint8_t val1 = (var_idx1 < c->num_pis) ? ((input_vector >> var_idx1) & 1) : gate_values[var_idx1 - c->num_pis];
        uint8_t final_val1 = val1 ^ is_inv1;

        // Evaluate second input literal
        uint32_t var_idx2 = g.input2 >> 1;
        uint8_t is_inv2 = g.input2 & 1;
        uint8_t val2 = (var_idx2 < c->num_pis) ? ((input_vector >> var_idx2) & 1) : gate_values[var_idx2 - c->num_pis];
        uint8_t final_val2 = val2 ^ is_inv2;

        gate_values[i] = final_val1 & final_val2;
    }

    // Evaluate the final circuit output literal
    uint32_t out_var_idx = c->output_literal >> 1;
    uint8_t out_is_inv = c->output_literal & 1;
    uint8_t out_val = (out_var_idx < c->num_pis) ? ((input_vector >> out_var_idx) & 1) : gate_values[out_var_idx - c->num_pis];
    
    return out_val ^ out_is_inv;
}

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <num_pi> <gates_a> <gates_b> <sim_vectors> <seed>\n", argv[0]);
        exit(1);
    }

    NUM_PRIMARY_INPUTS = atoi(argv[1]);
    NUM_GATES_CIRCUIT_A = atoi(argv[2]);
    NUM_GATES_CIRCUIT_B = atoi(argv[3]);
    NUM_SIMULATION_VECTORS = atoi(argv[4]);
    uint32_t seed = atoi(argv[5]);
    mt_seed(seed);

    if (NUM_PRIMARY_INPUTS > 64) {
        fprintf(stderr, "Error: num_primary_inputs cannot exceed 64 for this implementation.\n");
        exit(1);
    }

    g_circuit_a = create_random_circuit(NUM_PRIMARY_INPUTS, NUM_GATES_CIRCUIT_A);
    g_circuit_b = create_random_circuit(NUM_PRIMARY_INPUTS, NUM_GATES_CIRCUIT_B);
    if (!g_circuit_a || !g_circuit_b) {
        fprintf(stderr, "Error: Failed to allocate memory for circuits\n");
        exit(1);
    }

    g_simulation_vectors = (uint64_t*)malloc(NUM_SIMULATION_VECTORS * sizeof(uint64_t));
    if (!g_simulation_vectors) {
        fprintf(stderr, "Error: Failed to allocate memory for simulation vectors\n");
        exit(1);
    }

    for (uint32_t i = 0; i < NUM_SIMULATION_VECTORS; ++i) {
        uint64_t high = (uint64_t)mt_rand() << 32;
        uint64_t low = mt_rand();
        g_simulation_vectors[i] = high | low;
    }
}

int run_computation() {
    // Allocate scratchpads for simulation values
    uint8_t* gate_values_a = (uint8_t*)malloc(g_circuit_a->num_gates * sizeof(uint8_t));
    uint8_t* gate_values_b = (uint8_t*)malloc(g_circuit_b->num_gates * sizeof(uint8_t));
    if (!gate_values_a || !gate_values_b) {
        fprintf(stderr, "Error: Failed to allocate simulation scratchpads\n");
        free(gate_values_a);
        free(gate_values_b);
        exit(1);
    }

    int mismatch_count = 0;
    for (uint32_t i = 0; i < NUM_SIMULATION_VECTORS; ++i) {
        uint64_t vector = g_simulation_vectors[i];
        uint8_t out_a = evaluate_circuit(g_circuit_a, vector, gate_values_a);
        uint8_t out_b = evaluate_circuit(g_circuit_b, vector, gate_values_b);
        if (out_a != out_b) {
            mismatch_count++;
        }
    }

    free(gate_values_a);
    free(gate_values_b);
    return mismatch_count;
}

void cleanup() {
    if (g_circuit_a) {
        free(g_circuit_a->gates);
        free(g_circuit_a);
    }
    if (g_circuit_b) {
        free(g_circuit_b->gates);
        free(g_circuit_b);
    }
    free(g_simulation_vectors);
}

int main(int argc, char* argv[]) {
    struct timespec start, end;
    int final_result;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    final_result = run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%d\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
