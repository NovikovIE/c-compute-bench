#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

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

// Benchmark specific data structures

typedef struct {
    int driver_gate_index;
    float arrival_time; // Computed value: arrival time of the signal on this net
} Net;

typedef struct {
    int* input_net_indices;
    int num_inputs;
    int output_net_index;
    float intrinsic_delay; // Delay inherent to the gate
} Gate;

// Global variables for benchmark data
static int NUM_GATES;
static int NUM_NETS;
static Gate* gates = NULL;
static Net* nets = NULL;
static double final_result = 0.0;

#define MAX_INPUTS 8
#define CHANCE_OF_PRIMARY_INPUT 5 // in percent

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_gates> <num_nets> <seed>\n", argv[0]);
        exit(1);
    }

    NUM_GATES = atoi(argv[1]);
    NUM_NETS = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    if (NUM_GATES <= 0 || NUM_NETS <= 0) {
        fprintf(stderr, "Number of gates and nets must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    gates = (Gate*)malloc(NUM_GATES * sizeof(Gate));
    nets = (Net*)malloc(NUM_NETS * sizeof(Net));
    if (!gates || !nets) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(1);
    }

    // Initialize nets
    for (int i = 0; i < NUM_NETS; ++i) {
        nets[i].arrival_time = 0.0f;
        // Each net is driven by one gate, determined by modulo
        nets[i].driver_gate_index = i % NUM_GATES;
    }

    // Initialize gates with random properties and connections
    for (int i = 0; i < NUM_GATES; ++i) {
        // Each gate drives one net
        gates[i].output_net_index = i % NUM_NETS;
        gates[i].intrinsic_delay = (mt_rand() % 100) / 200.0f + 0.1f; // 0.1 to 0.6

        // Determine number of inputs. Some gates are primary inputs (0 inputs).
        if ((mt_rand() % 100) < CHANCE_OF_PRIMARY_INPUT) {
            gates[i].num_inputs = 0;
            gates[i].input_net_indices = NULL;
        } else {
            gates[i].num_inputs = 1 + (mt_rand() % MAX_INPUTS);
            gates[i].input_net_indices = (int*)malloc(gates[i].num_inputs * sizeof(int));
            if (!gates[i].input_net_indices) {
                fprintf(stderr, "Memory allocation for inputs failed.\n");
                exit(1);
            }
            for (int j = 0; j < gates[i].num_inputs; ++j) {
                // A gate can be connected to any net, potentially creating cycles.
                gates[i].input_net_indices[j] = mt_rand() % NUM_NETS;
            }
        }
    }
}

void run_computation() {
    int* in_degree = (int*)malloc(NUM_GATES * sizeof(int));
    int* queue = (int*)malloc(NUM_GATES * sizeof(int));
    int* sorted_gates = (int*)malloc(NUM_GATES * sizeof(int));

    if (!in_degree || !queue || !sorted_gates) {
        fprintf(stderr, "Failed to allocate temp memory for computation.\n");
        exit(1);
    }

    int queue_head = 0, queue_tail = 0;
    int sorted_count = 0;

    // 1. Initialize in-degrees and find primary inputs (gates with 0 inputs)
    for (int i = 0; i < NUM_GATES; ++i) {
        in_degree[i] = gates[i].num_inputs;
        if (in_degree[i] == 0) {
            queue[queue_tail++] = i;
        }
    }

    // 2. Perform Topological Sort (Kahn's algorithm) to find a valid evaluation order
    // This loop is the most computationally expensive part, as it finds the fanout
    // of each gate by iterating through all other gates.
    while (queue_head < queue_tail) {
        int u_idx = queue[queue_head++];
        sorted_gates[sorted_count++] = u_idx;

        int u_output_net = gates[u_idx].output_net_index;

        // Find all gates 'v' that are driven by gate 'u'
        for (int v_idx = 0; v_idx < NUM_GATES; ++v_idx) {
            for (int k = 0; k < gates[v_idx].num_inputs; ++k) {
                if (gates[v_idx].input_net_indices[k] == u_output_net) {
                    in_degree[v_idx]--;
                    if (in_degree[v_idx] == 0) {
                        queue[queue_tail++] = v_idx;
                    }
                    // Optimization: if a gate has multiple inputs from the same net,
                    // we only decrement once per gate, not per pin.
                    break;
                }
            }
        }
    }

    // 3. Propagate arrival times through the topologically sorted gates
    for (int i = 0; i < sorted_count; ++i) {
        int gate_idx = sorted_gates[i];

        float max_input_arrival_time = 0.0f;
        for (int j = 0; j < gates[gate_idx].num_inputs; ++j) {
            int input_net = gates[gate_idx].input_net_indices[j];
            if (nets[input_net].arrival_time > max_input_arrival_time) {
                max_input_arrival_time = nets[input_net].arrival_time;
            }
        }

        float output_arrival_time = max_input_arrival_time + gates[gate_idx].intrinsic_delay;
        int output_net = gates[gate_idx].output_net_index;
        // A net's arrival time is determined by the latest arriving signal from its driver
        if(output_arrival_time > nets[output_net].arrival_time) {
            nets[output_net].arrival_time = output_arrival_time;
        }
    }

    // 4. Calculate final result to prevent dead code elimination
    double total_arrival_time = 0.0;
    for (int i = 0; i < NUM_NETS; ++i) {
        total_arrival_time += nets[i].arrival_time;
    }
    final_result = total_arrival_time;

    free(in_degree);
    free(queue);
    free(sorted_gates);
}

void cleanup() {
    if (gates) {
        for (int i = 0; i < NUM_GATES; ++i) {
            if (gates[i].input_net_indices) {
                free(gates[i].input_net_indices);
            }
        }
        free(gates);
    }
    if (nets) {
        free(nets);
    }
    gates = NULL;
    nets = NULL;
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
    printf("%f\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
