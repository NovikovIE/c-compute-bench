/**
 * @file power_analysis_simulation.c
 * @brief A benchmark simulating power consumption in a digital circuit.
 * 
 * This program models a simplified digital circuit consisting of combinational gates
 * and flip-flops. It simulates the circuit's activity over a number of clock
 * cycles and calculates the total power consumed, which is comprised of dynamic
 * (switching) power and static (leakage) power.
 *
 * - Dynamic Power: Incurred when a gate or flip-flop changes state (toggles).
 *   It is modeled as being proportional to a pre-assigned switching capacitance.
 * - Static Power: Incurred by all components in every clock cycle, regardless
 *   of activity, modeling leakage current.
 *
 * The simulation is stochastic. The toggling of components is determined by a
 * 'toggle rate' probability, representing the overall activity factor of the circuit.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- Mersenne Twister (MT19937) PRNG --- VERBATIM AS PROVIDED ---
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

// Global structure to hold benchmark data and parameters
struct {
    // Parameters
    int num_gates;
    int num_flip_flops;
    int num_clock_cycles;
    float toggle_rate;

    // Circuit state and properties
    uint8_t *gate_states;
    uint8_t *flip_flop_states;
    float *gate_switching_capacitance;
    float *gate_leakage_power;

    // Result accumulator
    double total_power;
} g_data;

/**
 * @brief Parses arguments, allocates memory, and initializes circuit data.
 */
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <num_gates> <num_flip_flops> <num_clock_cycles> <toggle_rate_percentage> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_gates = atoi(argv[1]);
    g_data.num_flip_flops = atoi(argv[2]);
    g_data.num_clock_cycles = atoi(argv[3]);
    g_data.toggle_rate = atof(argv[4]) / 100.0f;
    uint32_t seed = (uint32_t)atoi(argv[5]);

    mt_seed(seed);

    // Allocate memory on the heap
    g_data.gate_states = (uint8_t *)malloc(g_data.num_gates * sizeof(uint8_t));
    g_data.gate_switching_capacitance = (float *)malloc(g_data.num_gates * sizeof(float));
    g_data.gate_leakage_power = (float *)malloc(g_data.num_gates * sizeof(float));
    g_data.flip_flop_states = (uint8_t *)malloc(g_data.num_flip_flops * sizeof(uint8_t));

    if (!g_data.gate_states || !g_data.gate_switching_capacitance || 
        !g_data.gate_leakage_power || !g_data.flip_flop_states) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize circuit properties and states randomly
    for (int i = 0; i < g_data.num_gates; ++i) {
        g_data.gate_states[i] = mt_rand() % 2;
        // Assign random capacitance (for dynamic power)
        g_data.gate_switching_capacitance[i] = (mt_rand() / (float)UINT32_MAX) * 0.9f + 0.1f;
        // Assign random leakage (for static power)
        g_data.gate_leakage_power[i] = (mt_rand() / (float)UINT32_MAX) * 0.009f + 0.001f;
    }

    for (int i = 0; i < g_data.num_flip_flops; ++i) {
        g_data.flip_flop_states[i] = mt_rand() % 2;
    }

    g_data.total_power = 0.0;
}

/**
 * @brief Runs the core power analysis simulation.
 */
void run_computation() {
    double total_dynamic_power = 0.0;
    double total_static_power = 0.0;

    for (int cycle = 0; cycle < g_data.num_clock_cycles; ++cycle) {
        // 1. Simulate combinational gate activity
        for (int i = 0; i < g_data.num_gates; ++i) {
            // Stochastically determine if a gate toggles based on the activity rate
            if ((mt_rand() / (float)UINT32_MAX) < g_data.toggle_rate) {
                g_data.gate_states[i] = 1 - g_data.gate_states[i]; // Toggle state
                total_dynamic_power += g_data.gate_switching_capacitance[i];
            }
            // Static power is consumed regardless of activity
            total_static_power += g_data.gate_leakage_power[i];
        }

        // 2. Simulate flip-flop activity
        // Here, their power contribution is simplified and rolled into the gate simulation,
        // but we still simulate their state changes as they drive the logic.
        for (int i = 0; i < g_data.num_flip_flops; ++i) {
             if ((mt_rand() / (float)UINT32_MAX) < g_data.toggle_rate) {
                g_data.flip_flop_states[i] = 1 - g_data.flip_flop_states[i];
            }
        }
    }
    g_data.total_power = total_dynamic_power + total_static_power;
}

/**
 * @brief Frees all memory allocated in setup_benchmark.
 */
void cleanup() {
    free(g_data.gate_states);
    free(g_data.gate_switching_capacitance);
    free(g_data.gate_leakage_power);
    free(g_data.flip_flop_states);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final accumulated power to stdout
    printf("%f\n", g_data.total_power);

    // Print the measured time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
