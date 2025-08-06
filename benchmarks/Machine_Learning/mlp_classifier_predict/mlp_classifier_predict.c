#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>
#include <string.h>

// --- Mersenne Twister (Do Not Modify) ---
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
// --- end Mersenne Twister ---

// Global parameters
int num_samples;
int num_features;
int num_hidden_layers;
int hidden_layer_size;
int num_classes;

// Global data structures for the MLP
float **input_samples; // [num_samples][num_features]
float ***weights;      // [num_hidden_layers + 1][rows][cols]
float **biases;        // [num_hidden_layers + 1][size]
long final_result;     // Accumulated result to prevent dead-code elimination

// Helper buffers for computation
float *activation_buffer_a;
float *activation_buffer_b;

// Helper to generate random float in [-1, 1]
float rand_float() {
    return ((float)mt_rand() / (float)UINT32_MAX) * 2.0f - 1.0f;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 7) {
        fprintf(stderr, "Usage: %s num_samples num_features num_hidden_layers hidden_layer_size num_classes seed\n", argv[0]);
        exit(1);
    }
    
    num_samples = atoi(argv[1]);
    num_features = atoi(argv[2]);
    num_hidden_layers = atoi(argv[3]);
    hidden_layer_size = atoi(argv[4]);
    num_classes = atoi(argv[5]);
    uint32_t seed = (uint32_t)atoi(argv[6]);
    
    mt_seed(seed);

    // Allocate input samples
    input_samples = (float **)malloc(num_samples * sizeof(float *));
    for (int i = 0; i < num_samples; ++i) {
        input_samples[i] = (float *)malloc(num_features * sizeof(float));
        for (int j = 0; j < num_features; ++j) {
            input_samples[i][j] = rand_float();
        }
    }

    int total_layers = num_hidden_layers + 1;
    weights = (float ***)malloc(total_layers * sizeof(float **));
    biases = (float **)malloc(total_layers * sizeof(float *));

    // Layer 0: input -> hidden 1
    int in_size = num_features;
    int out_size = (num_hidden_layers > 0) ? hidden_layer_size : num_classes;
    weights[0] = (float **)malloc(out_size * sizeof(float *));
    for (int i = 0; i < out_size; ++i) {
        weights[0][i] = (float *)malloc(in_size * sizeof(float));
        for (int j = 0; j < in_size; ++j) weights[0][i][j] = rand_float();
    }
    biases[0] = (float *)malloc(out_size * sizeof(float));
    for (int i = 0; i < out_size; ++i) biases[0][i] = rand_float();

    // Hidden layers: hidden -> hidden
    for (int l = 1; l < num_hidden_layers; ++l) {
        in_size = hidden_layer_size;
        out_size = hidden_layer_size;
        weights[l] = (float **)malloc(out_size * sizeof(float *));
        for (int i = 0; i < out_size; ++i) {
            weights[l][i] = (float *)malloc(in_size * sizeof(float));
            for (int j = 0; j < in_size; ++j) weights[l][i][j] = rand_float();
        }
        biases[l] = (float *)malloc(out_size * sizeof(float));
        for (int i = 0; i < out_size; ++i) biases[l][i] = rand_float();
    }

    // Output layer: last hidden -> output
    if (num_hidden_layers > 0) {
        in_size = hidden_layer_size;
        out_size = num_classes;
        weights[num_hidden_layers] = (float **)malloc(out_size * sizeof(float *));
        for (int i = 0; i < out_size; ++i) {
            weights[num_hidden_layers][i] = (float *)malloc(in_size * sizeof(float));
            for (int j = 0; j < in_size; ++j) weights[num_hidden_layers][i][j] = rand_float();
        }
        biases[num_hidden_layers] = (float *)malloc(out_size * sizeof(float));
        for (int i = 0; i < out_size; ++i) biases[num_hidden_layers][i] = rand_float();
    }

    // Allocate activation buffers for computation
    int max_layer_size = num_features > hidden_layer_size ? num_features : hidden_layer_size;
    max_layer_size = max_layer_size > num_classes ? max_layer_size : num_classes;
    if (max_layer_size == 0) max_layer_size = 1; // handle edge case of no layers

    activation_buffer_a = (float*)malloc(max_layer_size * sizeof(float));
    activation_buffer_b = (float*)malloc(max_layer_size * sizeof(float));
}

void run_computation() {
    long result_accumulator = 0;

    float *input_activations = activation_buffer_a;
    float *output_activations = activation_buffer_b;

    for (int s = 0; s < num_samples; ++s) {
        memcpy(input_activations, input_samples[s], num_features * sizeof(float));

        int total_layers = num_hidden_layers + 1;
        for (int l = 0; l < total_layers; ++l) {
            int in_size = (l == 0) ? num_features : hidden_layer_size;
            int out_size = (l == num_hidden_layers) ? num_classes : hidden_layer_size;

            for (int i = 0; i < out_size; ++i) {
                float sum = biases[l][i];
                for (int j = 0; j < in_size; ++j) {
                    sum += weights[l][i][j] * input_activations[j];
                }
                output_activations[i] = sum;
            }
            
            if (l < num_hidden_layers) { 
                for (int i = 0; i < out_size; ++i) {
                   output_activations[i] = tanhf(output_activations[i]);
                }
            }

            float *temp = input_activations;
            input_activations = output_activations;
            output_activations = temp;
        }
        
        float *final_activations = input_activations;
        float max_val = -INFINITY;
        int prediction = 0;
        for (int i = 0; i < num_classes; ++i) {
            if (final_activations[i] > max_val) {
                max_val = final_activations[i];
                prediction = i;
            }
        }
        result_accumulator += prediction;
    }
    final_result = result_accumulator;
}

void cleanup() {
    for (int i = 0; i < num_samples; ++i) {
        free(input_samples[i]);
    }
    free(input_samples);

    int total_layers = num_hidden_layers + 1;
    for (int l = 0; l < total_layers; ++l) {
        int out_size = (l == num_hidden_layers) ? num_classes : hidden_layer_size;
        // For layer 0, if there are no hidden layers, out_size is num_classes
        if (l == 0 && num_hidden_layers == 0) out_size = num_classes;
        
        for (int i = 0; i < out_size; ++i) {
            free(weights[l][i]);
        }
        free(weights[l]);
        free(biases[l]);
    }
    free(weights);
    free(biases);

    free(activation_buffer_a);
    free(activation_buffer_b);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%ld\n", final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
