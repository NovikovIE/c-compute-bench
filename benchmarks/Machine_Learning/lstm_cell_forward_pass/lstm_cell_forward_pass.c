#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- Mersenne Twister (MT19937) Generator ---
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
// --- End of Mersenne Twister ---

// Benchmark parameters
int batch_size;
int input_features;
int output_features;
int combined_features;

// Data arrays
float *x_t;           // Input
float *h_prev;        // Previous hidden state
float *c_prev;        // Previous cell state
float *W_all;         // Combined weights for all 4 gates
float *b_all;         // Combined biases for all 4 gates
float *h_next;        // Next hidden state (output)
float *c_next;        // Next cell state (output)

// Intermediate buffers for computation
float *combined_input;
float *pre_activations;

// Result accumulator
double final_result;

// Helper to generate a random float between -1.0 and 1.0
float gen_random_float() {
    return ((float)mt_rand() / (float)UINT32_MAX) * 2.0f - 1.0f;
}

// Sigmoid activation function
float sigmoidf(float x) {
    return 1.0f / (1.0f + expf(-x));
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <batch_size> <input_features> <output_features> <seed>\n", argv[0]);
        exit(1);
    }

    batch_size = atoi(argv[1]);
    input_features = atoi(argv[2]);
    output_features = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    if (batch_size <= 0 || input_features <= 0 || output_features <= 0) {
        fprintf(stderr, "FATAL: all dimensions must be positive integers.\n");
        exit(1);
    }
    
    mt_seed(seed);

    combined_features = input_features + output_features;

    // Allocate memory for data and parameters
    x_t = (float*)malloc(batch_size * input_features * sizeof(float));
    h_prev = (float*)malloc(batch_size * output_features * sizeof(float));
    c_prev = (float*)malloc(batch_size * output_features * sizeof(float));
    W_all = (float*)malloc(combined_features * (4 * output_features) * sizeof(float));
    b_all = (float*)malloc(4 * output_features * sizeof(float));
    h_next = (float*)malloc(batch_size * output_features * sizeof(float));
    c_next = (float*)malloc(batch_size * output_features * sizeof(float));
    combined_input = (float*)malloc(batch_size * combined_features * sizeof(float));
    pre_activations = (float*)malloc(batch_size * (4 * output_features) * sizeof(float));

    if (!x_t || !h_prev || !c_prev || !W_all || !b_all || !h_next || !c_next || !combined_input || !pre_activations) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize data with random values
    for (int i = 0; i < batch_size * input_features; i++) x_t[i] = gen_random_float();
    for (int i = 0; i < batch_size * output_features; i++) h_prev[i] = gen_random_float();
    for (int i = 0; i < batch_size * output_features; i++) c_prev[i] = gen_random_float();
    for (int i = 0; i < combined_features * (4 * output_features); i++) W_all[i] = gen_random_float() * 0.1f;
    for (int i = 0; i < 4 * output_features; i++) b_all[i] = gen_random_float() * 0.1f;
}

void run_computation() {
    // 1. Prepare combined input: concatenate x_t and h_prev for each batch item
    for (int b = 0; b < batch_size; ++b) {
        // Copy x_t part
        for (int i = 0; i < input_features; ++i) {
            combined_input[b * combined_features + i] = x_t[b * input_features + i];
        }
        // Copy h_prev part
        for (int i = 0; i < output_features; ++i) {
            combined_input[b * combined_features + input_features + i] = h_prev[b * output_features + i];
        }
    }

    // 2. Main matrix multiplication: pre_activations = combined_input @ W_all
    //  (batch_size, combined_features) @ (combined_features, 4 * output_features) -> (batch_size, 4 * output_features)
    for (int b = 0; b < batch_size; ++b) {
        for (int j = 0; j < (4 * output_features); ++j) {
            float sum = 0.0f;
            for (int k = 0; k < combined_features; ++k) {
                sum += combined_input[b * combined_features + k] * W_all[k * (4 * output_features) + j];
            }
            pre_activations[b * (4 * output_features) + j] = sum;
        }
    }

    // 3. Element-wise operations: biases, activations, and state updates
    for (int b = 0; b < batch_size; ++b) {
        for (int j = 0; j < output_features; ++j) {
            float* pre_act_base = &pre_activations[b * (4 * output_features)];
            
            float pre_i = pre_act_base[j]                             + b_all[j];
            float pre_f = pre_act_base[j + output_features]           + b_all[j + output_features];
            float pre_g = pre_act_base[j + 2 * output_features]       + b_all[j + 2 * output_features];
            float pre_o = pre_act_base[j + 3 * output_features]       + b_all[j + 3 * output_features];

            float i_t = sigmoidf(pre_i);
            float f_t = sigmoidf(pre_f);
            float g_t = tanhf(pre_g);
            float o_t = sigmoidf(pre_o);

            int state_idx = b * output_features + j;
            c_next[state_idx] = (f_t * c_prev[state_idx]) + (i_t * g_t);
            h_next[state_idx] = o_t * tanhf(c_next[state_idx]);
        }
    }

    // 4. Accumulate result to prevent dead code elimination
    final_result = 0.0;
    for (int i = 0; i < batch_size * output_features; ++i) {
        final_result += h_next[i];
    }
}

void cleanup() {
    free(x_t);
    free(h_prev);
    free(c_prev);
    free(W_all);
    free(b_all);
    free(h_next);
    free(c_next);
    free(combined_input);
    free(pre_activations);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%f\n", final_result);

    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
