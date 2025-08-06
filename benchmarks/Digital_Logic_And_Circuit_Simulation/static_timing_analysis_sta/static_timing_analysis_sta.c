#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- MERSENNE TWISTER (Verbatim) ---
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

// --- BENCHMARK DATA STRUCTURES ---

typedef struct {
    double arrival_time;
    double gate_delay;
    int fan_out_count;
    int* fan_out_indices;
} Gate;

// Global structure to hold benchmark data
struct {
    int num_gates;
    Gate* gates;
    int* fanout_buffer;
    long long final_result;
} benchmark_data;

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_gates> <edge_factor> <seed>\n", argv[0]);
        exit(1);
    }

    benchmark_data.num_gates = atoi(argv[1]);
    int edge_factor = atoi(argv[2]); // Average fan-out per gate
    uint32_t seed = atoi(argv[3]);

    if (benchmark_data.num_gates <= 0 || edge_factor < 0) {
        fprintf(stderr, "ERROR: num_gates must be positive, edge_factor non-negative.\n");
        exit(1);
    }

    mt_seed(seed);

    benchmark_data.gates = (Gate*)malloc(benchmark_data.num_gates * sizeof(Gate));
    if (!benchmark_data.gates) {
        fprintf(stderr, "ERROR: Failed to allocate memory for gates.\n");
        exit(1);
    }
    
    long long total_edges = (long long)benchmark_data.num_gates * edge_factor;
    if (total_edges == 0) total_edges = 1;

    benchmark_data.fanout_buffer = (int*)malloc(total_edges * sizeof(int));
    if (!benchmark_data.fanout_buffer) {
        fprintf(stderr, "ERROR: Failed to allocate memory for fanout buffer.\n");
        free(benchmark_data.gates);
        exit(1);
    }

    int* current_fanout_ptr = benchmark_data.fanout_buffer;
    
    for (int i = 0; i < benchmark_data.num_gates; i++) {
        Gate* current_gate = &benchmark_data.gates[i];
        
        current_gate->arrival_time = 0.0;
        current_gate->gate_delay = (double)(mt_rand() % 20) + 1.0; // Delay of 1 to 20 units

        int fan_out_count = 0;
        if (i < benchmark_data.num_gates - 1 && edge_factor > 0) {
            fan_out_count = edge_factor;
        }

        if ((current_fanout_ptr - benchmark_data.fanout_buffer + fan_out_count) > total_edges) {
            fan_out_count = total_edges - (current_fanout_ptr - benchmark_data.fanout_buffer);
        }

        current_gate->fan_out_count = fan_out_count;
        current_gate->fan_out_indices = current_fanout_ptr;
        
        if (fan_out_count > 0) {
            int range = benchmark_data.num_gates - (i + 1);
            for (int j = 0; j < fan_out_count; j++) {
                int target_idx = (mt_rand() % range) + (i + 1);
                current_gate->fan_out_indices[j] = target_idx;
            }
        }
        current_fanout_ptr += fan_out_count;
    }

    benchmark_data.final_result = 0;
}

void run_computation() {
    for (int i = 0; i < benchmark_data.num_gates; i++) {
        Gate* current_gate = &benchmark_data.gates[i];

        double output_arrival_time = current_gate->arrival_time + current_gate->gate_delay;
        
        for (int j = 0; j < current_gate->fan_out_count; j++) {
            Gate* next_gate = &benchmark_data.gates[current_gate->fan_out_indices[j]];
            
            if (next_gate->arrival_time < output_arrival_time) {
                next_gate->arrival_time = output_arrival_time;
            }
        }
    }

    double total_arrival_time = 0.0;
    for (int i = 0; i < benchmark_data.num_gates; i++) {
        total_arrival_time += benchmark_data.gates[i].arrival_time;
    }
    benchmark_data.final_result = (long long)total_arrival_time;
}

void cleanup() {
    free(benchmark_data.gates);
    free(benchmark_data.fanout_buffer);
    benchmark_data.gates = NULL;
    benchmark_data.fanout_buffer = NULL;
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
    printf("%lld\n", benchmark_data.final_result);
    
    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}
