#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- START MERSENNE TWISTER (MT19937) --- Do Not Modify ---
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

// --- BENCHMARK DATA & PARAMETERS ---
typedef struct {
    // Parameters
    int num_samples;
    int num_features;
    int num_hidden_layers;
    int hidden_layer_size;
    int num_output_neurons;
    float learning_rate;

    // Data
    float** inputs;
    float** targets;

    // Network architecture
    int num_total_layers;         // num_hidden_layers + 2 (input, hidden(s), output)
    int* layer_sizes;

    // Network state
    float*** weights;             // Array of weight matrices
    float** biases;               // Array of bias vectors
    float** activations;          // Array of activation vectors (neuron outputs)
    float** deltas;               // Array of error gradient vectors

    // Final result accumulator
    double final_result;
} BenchmarkData;

static BenchmarkData g_data;

// --- HELPER FUNCTIONS ---
float rand_float() {
    return ((float)mt_rand() / (float)UINT32_MAX) * 2.0f - 1.0f;
}

float sigmoid(float x) {
    return 1.0f / (1.0f + expf(-x));
}

// --- BENCHMARK CORE FUNCTIONS ---

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 8) {
        fprintf(stderr, "Usage: %s num_samples num_features num_hidden_layers hidden_layer_size num_output_neurons learning_rate seed\n", argv[0]);
        exit(1);
    }

    g_data.num_samples = atoi(argv[1]);
    g_data.num_features = atoi(argv[2]);
    g_data.num_hidden_layers = atoi(argv[3]);
    g_data.hidden_layer_size = atoi(argv[4]);
    g_data.num_output_neurons = atoi(argv[5]);
    g_data.learning_rate = atof(argv[6]);
    uint32_t seed = (uint32_t)atoi(argv[7]);

    mt_seed(seed);

    // --- Define network architecture ---
    g_data.num_total_layers = g_data.num_hidden_layers + 2;
    g_data.layer_sizes = (int*)malloc(g_data.num_total_layers * sizeof(int));
    if (!g_data.layer_sizes) { perror("malloc failed"); exit(1); }

    g_data.layer_sizes[0] = g_data.num_features;
    for (int i = 1; i <= g_data.num_hidden_layers; ++i) {
        g_data.layer_sizes[i] = g_data.hidden_layer_size;
    }
    g_data.layer_sizes[g_data.num_total_layers - 1] = g_data.num_output_neurons;

    // --- Allocate and initialize input/target data ---
    g_data.inputs = (float**)malloc(g_data.num_samples * sizeof(float*));
    g_data.targets = (float**)malloc(g_data.num_samples * sizeof(float*));
    if (!g_data.inputs || !g_data.targets) { perror("malloc failed"); exit(1); }

    for (int i = 0; i < g_data.num_samples; ++i) {
        g_data.inputs[i] = (float*)malloc(g_data.num_features * sizeof(float));
        g_data.targets[i] = (float*)malloc(g_data.num_output_neurons * sizeof(float));
        if (!g_data.inputs[i] || !g_data.targets[i]) { perror("malloc failed"); exit(1); }
        for (int j = 0; j < g_data.num_features; ++j) g_data.inputs[i][j] = rand_float();
        for (int j = 0; j < g_data.num_output_neurons; ++j) g_data.targets[i][j] = rand_float();
    }

    // --- Allocate network structures ---
    int num_comp_layers = g_data.num_total_layers - 1;
    g_data.weights = (float***)malloc(num_comp_layers * sizeof(float**));
    g_data.biases = (float**)malloc(num_comp_layers * sizeof(float*));
    g_data.deltas = (float**)malloc(num_comp_layers * sizeof(float*));
    g_data.activations = (float**)malloc(g_data.num_total_layers * sizeof(float*));
    if (!g_data.weights || !g_data.biases || !g_data.deltas || !g_data.activations) { perror("malloc failed"); exit(1); }

    for (int i = 0; i < g_data.num_total_layers; ++i) {
        g_data.activations[i] = (float*)calloc(g_data.layer_sizes[i], sizeof(float));
        if (!g_data.activations[i]) { perror("calloc failed"); exit(1); }
    }

    for (int i = 0; i < num_comp_layers; ++i) {
        int in_size = g_data.layer_sizes[i];
        int out_size = g_data.layer_sizes[i + 1];
        g_data.weights[i] = (float**)malloc(out_size * sizeof(float*));
        g_data.biases[i] = (float*)malloc(out_size * sizeof(float));
        g_data.deltas[i] = (float*)malloc(out_size * sizeof(float));
        if (!g_data.weights[i] || !g_data.biases[i] || !g_data.deltas[i]) { perror("malloc failed"); exit(1); }

        for (int j = 0; j < out_size; ++j) {
            g_data.weights[i][j] = (float*)malloc(in_size * sizeof(float));
            if (!g_data.weights[i][j]) { perror("malloc failed"); exit(1); }
            for (int k = 0; k < in_size; ++k) {
                g_data.weights[i][j][k] = rand_float();
            }
            g_data.biases[i][j] = rand_float();
        }
    }
    g_data.final_result = 0.0;
}

void run_computation() {
    int num_comp_layers = g_data.num_total_layers - 1;

    for (int s = 0; s < g_data.num_samples; ++s) {
        // Set input layer activations to the current sample's features
        for (int i = 0; i < g_data.num_features; ++i) {
            g_data.activations[0][i] = g_data.inputs[s][i];
        }

        // --- 1. Feedforward pass ---
        for (int l = 0; l < num_comp_layers; ++l) {
            int in_size = g_data.layer_sizes[l];
            int out_size = g_data.layer_sizes[l + 1];
            for (int j = 0; j < out_size; ++j) {
                float z = g_data.biases[l][j];
                for (int k = 0; k < in_size; ++k) {
                    z += g_data.weights[l][j][k] * g_data.activations[l][k];
                }
                g_data.activations[l + 1][j] = sigmoid(z);
            }
        }

        // --- 2. Backpropagation pass ---
        // Calculate output layer deltas
        int out_layer_idx = num_comp_layers - 1;
        int out_size = g_data.layer_sizes[out_layer_idx + 1];
        for (int j = 0; j < out_size; ++j) {
            float activation = g_data.activations[out_layer_idx + 1][j];
            float error = g_data.targets[s][j] - activation;
            g_data.deltas[out_layer_idx][j] = error * activation * (1.0f - activation); // sigmoid derivative
        }

        // Propagate deltas to hidden layers
        for (int l = num_comp_layers - 2; l >= 0; --l) {
            int curr_size = g_data.layer_sizes[l + 1];
            int next_size = g_data.layer_sizes[l + 2];
            for (int j = 0; j < curr_size; ++j) {
                float error = 0.0f;
                for (int k = 0; k < next_size; ++k) {
                     error += g_data.weights[l + 1][k][j] * g_data.deltas[l + 1][k];
                }
                float activation = g_data.activations[l + 1][j];
                g_data.deltas[l][j] = error * activation * (1.0f - activation); // sigmoid derivative
            }
        }

        // --- 3. Update weights and biases ---
        for (int l = 0; l < num_comp_layers; ++l) {
            int in_size = g_data.layer_sizes[l];
            int out_size = g_data.layer_sizes[l + 1];
            for (int j = 0; j < out_size; ++j) {
                for (int k = 0; k < in_size; ++k) {
                    g_data.weights[l][j][k] += g_data.learning_rate * g_data.deltas[l][j] * g_data.activations[l][k];
                }
                g_data.biases[l][j] += g_data.learning_rate * g_data.deltas[l][j];
            }
        }
    }

    // --- Calculate final result to prevent dead code elimination ---
    double sum = 0.0;
    int final_layer_idx = g_data.num_total_layers - 1;
    int final_layer_size = g_data.layer_sizes[final_layer_idx];
    for (int i = 0; i < final_layer_size; ++i) {
        sum += g_data.activations[final_layer_idx][i];
    }
    g_data.final_result = sum;
}

void cleanup() {
    // Free input/target data
    for (int i = 0; i < g_data.num_samples; ++i) {
        free(g_data.inputs[i]);
        free(g_data.targets[i]);
    }
    free(g_data.inputs);
    free(g_data.targets);

    // Free network structures
    int num_comp_layers = g_data.num_total_layers - 1;
    for (int i = 0; i < num_comp_layers; ++i) {
        int out_size = g_data.layer_sizes[i + 1];
        for (int j = 0; j < out_size; ++j) {
            free(g_data.weights[i][j]);
        }
        free(g_data.weights[i]);
        free(g_data.biases[i]);
        free(g_data.deltas[i]);
    }
    free(g_data.weights);
    free(g_data.biases);
    free(g_data.deltas);

    for (int i = 0; i < g_data.num_total_layers; ++i) {
        free(g_data.activations[i]);
    }
    free(g_data.activations);

    free(g_data.layer_sizes);
}

int main(int argc, char* argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%f\n", g_data.final_result);
    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}
