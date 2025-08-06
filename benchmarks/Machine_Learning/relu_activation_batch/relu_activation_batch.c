/**
 * @file relu_activation_batch.c
 * @brief Benchmark for the ReLU (Rectified Linear Unit) activation function.
 *
 * This program benchmarks the performance of applying the ReLU activation
 * function to a batch of neuron outputs. The ReLU function is a fundamental
 * component in many neural networks, defined as f(x) = max(0, x).
 *
 * The benchmark initializes a large matrix of floating-point numbers representing
 * the pre-activation outputs for a batch of data instances across many neurons.
 * The core computation then applies the ReLU function to each element.
 *
 * The program is structured to separate data setup, computation, and cleanup.
 * Timing is performed only on the core computation part to ensure accuracy.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- Mersenne Twister (MT19937) Generator (Verbatim) ---
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

// --- Global Benchmark Data ---

// Parameters
static int g_num_neurons;
static int g_batch_size;

// Data arrays (flattened for contiguous memory access)
static float *g_inputs;
static float *g_outputs;

// Result accumulator to prevent dead code elimination
static double g_result_sum;

/**
 * @brief Generates a random float between -1.0 and 1.0.
 * Uses the MT19937 generator.
 */
float random_float() {
    // Scale the 32-bit unsigned integer to a float in [0, 1], then shift to [-1, 1]
    return 2.0f * ((float)mt_rand() / (float)UINT32_MAX) - 1.0f;
}

/**
 * @brief Sets up the benchmark by parsing arguments, allocating memory,
 *        and initializing data with random values.
 */
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_neurons> <batch_size> <seed>\n", argv[0]);
        exit(1);
    }

    g_num_neurons = atoi(argv[1]);
    g_batch_size = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (g_num_neurons <= 0 || g_batch_size <= 0) {
        fprintf(stderr, "FATAL: num_neurons and batch_size must be positive integers.\n");
        exit(1);
    }

    mt_seed(seed);

    size_t total_elements = (size_t)g_batch_size * g_num_neurons;
    g_inputs = (float *)malloc(total_elements * sizeof(float));
    g_outputs = (float *)malloc(total_elements * sizeof(float));

    if (g_inputs == NULL || g_outputs == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed for %zu elements.\n", total_elements);
        exit(1);
    }

    for (size_t i = 0; i < total_elements; ++i) {
        g_inputs[i] = random_float();
    }

    g_result_sum = 0.0;
}

/**
 * @brief Runs the core computation: applying ReLU to the input batch.
 */
void run_computation() {
    size_t total_elements = (size_t)g_batch_size * g_num_neurons;
    double local_sum = 0.0; // Use local accumulator for performance

    for (size_t i = 0; i < total_elements; ++i) {
        float input_val = g_inputs[i];
        // ReLU function: f(x) = max(0, x)
        float output_val = (input_val > 0.0f) ? input_val : 0.0f;
        g_outputs[i] = output_val;
        local_sum += output_val;
    }

    g_result_sum = local_sum;
}

/**
 * @brief Frees all memory allocated during setup.
 */
void cleanup() {
    free(g_inputs);
    free(g_outputs);
}

/**
 * @brief Main function to orchestrate the benchmark.
 */
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final accumulated result to stdout to prevent optimizer elimination
    printf("%f\n", g_result_sum);

    // Print the timing information to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
