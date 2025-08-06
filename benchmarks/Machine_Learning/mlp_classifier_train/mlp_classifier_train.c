#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>
#include <string.h>

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

// --- Global Data Structures ---
// Parameters
int num_samples;
int num_features;
int num_hidden_layers;
int hidden_layer_size;
int num_epochs;
const int num_classes = 10; // Fixed for simplicity
const float learning_rate = 0.01f;

// Data
float* X_data; // [num_samples x num_features]
int* y_data;   // [num_samples]

// MLP Model Structure
int* layer_sizes;
float** weights;
float** biases;
float** activations;
float** deltas;
float** pre_activations; // z values before activation function

// Result accumulator
float final_loss;

// --- Helper Functions ---
float rand_float() {
    return (mt_rand() / (float)UINT32_MAX) * 2.0f - 1.0f;
}

float relu(float x) {
    return x > 0 ? x : 0;
}

float relu_derivative(float x) {
    return x > 0 ? 1.0f : 0.0f;
}

// In-place softmax for the output layer's pre-activations
void softmax(float* arr, int size) {
    float max_val = arr[0];
    for (int i = 1; i < size; i++) {
        if (arr[i] > max_val) max_val = arr[i];
    }
    float sum_exp = 0.0f;
    for (int i = 0; i < size; i++) {
        arr[i] = expf(arr[i] - max_val); // Numerically stable
        sum_exp += arr[i];
    }
    for (int i = 0; i < size; i++) {
        arr[i] /= sum_exp;
    }
}

// --- Benchmark Functions ---
void setup_benchmark(int argc, char* argv[]) {
    if (argc != 7) {
        fprintf(stderr, "Usage: %s num_samples num_features num_hidden_layers hidden_layer_size num_epochs seed\n", argv[0]);
        exit(1);
    }

    num_samples = atoi(argv[1]);
    num_features = atoi(argv[2]);
    num_hidden_layers = atoi(argv[3]);
    hidden_layer_size = atoi(argv[4]);
    num_epochs = atoi(argv[5]);
    mt_seed((uint32_t)atoi(argv[6]));

    X_data = (float*)malloc(num_samples * num_features * sizeof(float));
    y_data = (int*)malloc(num_samples * sizeof(int));
    for (int i = 0; i < num_samples * num_features; i++) {
        X_data[i] = rand_float() * 0.5f;
    }
    for (int i = 0; i < num_samples; i++) {
        y_data[i] = mt_rand() % num_classes;
    }

    int num_total_layers = num_hidden_layers + 2;
    layer_sizes = (int*)malloc(num_total_layers * sizeof(int));
    layer_sizes[0] = num_features;
    for (int i = 1; i <= num_hidden_layers; i++) {
        layer_sizes[i] = hidden_layer_size;
    }
    layer_sizes[num_total_layers - 1] = num_classes;

    int num_weight_sets = num_total_layers - 1;
    weights = (float**)malloc(num_weight_sets * sizeof(float*));
    biases = (float**)malloc(num_weight_sets * sizeof(float*));
    deltas = (float**)malloc(num_weight_sets * sizeof(float*));
    pre_activations = (float**)malloc(num_weight_sets * sizeof(float*));
    activations = (float**)malloc(num_total_layers * sizeof(float*));

    for (int i = 0; i < num_weight_sets; i++) {
        int in_dim = layer_sizes[i];
        int out_dim = layer_sizes[i+1];
        weights[i] = (float*)malloc(in_dim * out_dim * sizeof(float));
        biases[i] = (float*)malloc(out_dim * sizeof(float));
        deltas[i] = (float*)malloc(out_dim * sizeof(float));
        pre_activations[i] = (float*)malloc(out_dim * sizeof(float));
        
        float limit = sqrtf(6.0f / (in_dim + out_dim));
        for (int j = 0; j < in_dim * out_dim; j++) {
            weights[i][j] = rand_float() * limit;
        }
        for (int j = 0; j < out_dim; j++) {
            biases[i][j] = 0.0f;
        }
    }

    for (int i = 0; i < num_total_layers; i++) {
        activations[i] = (float*)malloc(layer_sizes[i] * sizeof(float));
    }
}

void run_computation() {
    int num_total_layers = num_hidden_layers + 2;
    int num_weight_sets = num_total_layers - 1;
    float epoch_loss = 0.0f;

    for (int epoch = 0; epoch < num_epochs; epoch++) {
        epoch_loss = 0.0f;
        for (int s = 0; s < num_samples; s++) {
            // Forward Pass
            memcpy(activations[0], &X_data[s * num_features], num_features * sizeof(float));
            for (int l = 0; l < num_weight_sets; l++) {
                int in_dim = layer_sizes[l];
                int out_dim = layer_sizes[l + 1];
                for (int j = 0; j < out_dim; j++) {
                    float z = biases[l][j];
                    for (int i = 0; i < in_dim; i++) {
                        z += activations[l][i] * weights[l][i * out_dim + j];
                    }
                    pre_activations[l][j] = z;
                }
                if (l < num_weight_sets - 1) { // Hidden layers
                    for (int j = 0; j < out_dim; j++) {
                        activations[l + 1][j] = relu(pre_activations[l][j]);
                    }
                } else { // Output layer
                    memcpy(activations[l + 1], pre_activations[l], out_dim * sizeof(float));
                    softmax(activations[l + 1], out_dim);
                }
            }
            
            epoch_loss += -logf(activations[num_total_layers - 1][y_data[s]] + 1e-9f);

            // Backward Pass
            int output_layer_idx = num_weight_sets - 1;
            int out_dim = layer_sizes[output_layer_idx + 1];
            for (int i = 0; i < out_dim; i++) {
                deltas[output_layer_idx][i] = activations[output_layer_idx + 1][i];
            }
            deltas[output_layer_idx][y_data[s]] -= 1.0f;

            for (int l = num_weight_sets - 2; l >= 0; l--) {
                int current_dim = layer_sizes[l + 1];
                int next_dim = layer_sizes[l + 2];
                for (int j = 0; j < current_dim; j++) {
                    float error = 0.0f;
                    for (int k = 0; k < next_dim; k++) {
                        error += deltas[l + 1][k] * weights[l + 1][j * next_dim + k];
                    }
                    deltas[l][j] = error * relu_derivative(pre_activations[l][j]);
                }
            }

            // Update Weights and Biases
            for (int l = 0; l < num_weight_sets; l++) {
                int in_dim = layer_sizes[l];
                int out_dim = layer_sizes[l + 1];
                for (int j = 0; j < out_dim; j++) {
                    biases[l][j] -= learning_rate * deltas[l][j];
                    for (int i = 0; i < in_dim; i++) {
                        weights[l][i * out_dim + j] -= learning_rate * activations[l][i] * deltas[l][j];
                    }
                }
            }
        }
    }
    final_loss = epoch_loss / num_samples;
}

void cleanup() {
    int num_total_layers = num_hidden_layers + 2;
    int num_weight_sets = num_total_layers - 1;

    free(X_data);
    free(y_data);

    for (int i = 0; i < num_weight_sets; i++) {
        free(weights[i]);
        free(biases[i]);
        free(deltas[i]);
        free(pre_activations[i]);
    }
    free(weights);
    free(biases);
    free(deltas);
    free(pre_activations);

    for (int i = 0; i < num_total_layers; i++) {
        free(activations[i]);
    }
    free(activations);
    free(layer_sizes);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%f\n", final_loss);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
