#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- Mersenne Twister (MT19937) Generator (Do Not Modify) ---
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
// --- End of MT19937 --- 

// --- Benchmark Data Structures ---
typedef enum {
    G_IN, G_NOT, G_AND, G_OR, G_XOR, G_NAND, G_NOR, NUM_GATE_TYPES
} GateType;

typedef struct {
    GateType type;
    int num_inputs;
    int* inputs; // Array of indices to other gates
    int logic_level;
    double area_cost;
} Gate;

// --- Global Pointers & Parameters ---
static int num_initial_gates;
static int optimization_effort_level;

static Gate* circuit;
static long long final_result;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_initial_gates> <optimization_effort_level> <seed>\n", argv[0]);
        exit(1);
    }

    num_initial_gates = atoi(argv[1]);
    optimization_effort_level = atoi(argv[2]);
    uint32_t seed = (uint32_t)strtoul(argv[3], NULL, 10);

    if (num_initial_gates <= 0 || optimization_effort_level <= 0) {
        fprintf(stderr, "Parameters must be positive integers.\n");
        exit(1);
    }

    mt_seed(seed);

    circuit = (Gate*)malloc(num_initial_gates * sizeof(Gate));
    if (!circuit) {
        fprintf(stderr, "Failed to allocate memory for circuit\n");
        exit(1);
    }

    // Establish a base of primary inputs to prevent invalid circuits
    int num_primary_inputs = (num_initial_gates > 20) ? 20 : num_initial_gates;

    for (int i = 0; i < num_initial_gates; ++i) {
        circuit[i].logic_level = 0;
        
        if (i < num_primary_inputs) {
            circuit[i].type = G_IN;
            circuit[i].num_inputs = 0;
            circuit[i].inputs = NULL;
        } else {
            circuit[i].type = (GateType)((mt_rand() % (NUM_GATE_TYPES - 1)) + 1); // Avoid generating new G_IN gates
            circuit[i].num_inputs = (circuit[i].type == G_NOT) ? 1 : 2 + (mt_rand() % 2); // 1 input for NOT, 2-3 for others
            circuit[i].inputs = (int*)malloc(circuit[i].num_inputs * sizeof(int));
            if (!circuit[i].inputs) {
                fprintf(stderr, "Failed to allocate memory for gate inputs\n");
                exit(1);
            }
            
            int max_level = 0;
            for (int j = 0; j < circuit[i].num_inputs; ++j) {
                int input_idx = mt_rand() % i; // Guarantees a Directed Acyclic Graph (DAG)
                circuit[i].inputs[j] = input_idx;
                if (circuit[input_idx].logic_level > max_level) {
                    max_level = circuit[input_idx].logic_level;
                }
            }
            circuit[i].logic_level = max_level + 1;
        }

        // Assign a base cost. Cost is a simple function of inputs and type.
        circuit[i].area_cost = circuit[i].num_inputs * 1.5 + circuit[i].type * 0.5;
    }
}

void run_computation() {
    double total_area_saved = 0.0;
    const int search_attempts = 5; // How many alternative implementations to check for each gate

    // The optimization effort level controls the number of optimization passes
    for (int pass = 0; pass < optimization_effort_level; ++pass) {
        // In each pass, attempt to optimize every gate
        for (int i = 0; i < num_initial_gates; ++i) {
            Gate* g = &circuit[i];

            if (g->type == G_IN) continue;

            double initial_cost = g->area_cost;
            double best_new_cost = initial_cost;

            // Simulate a local search for a better (cheaper) gate implementation.
            // This is the core computational loop of the benchmark.
            for (int k = 0; k < search_attempts; ++k) {
                GateType new_type = (GateType)((mt_rand() % (NUM_GATE_TYPES - 1)) + 1);
                
                // New cost is a dummy calculation but represents a realistic workload
                double hypothetical_cost = g->num_inputs * 1.5 + new_type * 0.5;
                hypothetical_cost -= (double)(mt_rand() % 1000) / 5000.0; // Add noise for variability

                if (hypothetical_cost < best_new_cost) {
                    best_new_cost = hypothetical_cost;
                }
            }
            
            if (best_new_cost < initial_cost) {
                total_area_saved += (initial_cost - best_new_cost);
                g->area_cost = best_new_cost; // "Apply" the optimization by updating the cost
            }
        }
    }
    
    // Accumulate a final checksum to prevent dead code elimination.
    // Sum the final costs and the total savings.
    double final_checksum = 0.0;
    for (int i = 0; i < num_initial_gates; i++) {
        final_checksum += circuit[i].area_cost;
    }

    final_result = (long long)(total_area_saved + final_checksum);
}

void cleanup() {
    for (int i = 0; i < num_initial_gates; ++i) {
        if (circuit[i].inputs != NULL) {
            free(circuit[i].inputs);
        }
    }
    free(circuit);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%lld\n", final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
