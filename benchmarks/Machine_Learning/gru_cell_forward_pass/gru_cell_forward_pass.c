#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- START MERSENNE TWISTER (DO NOT MODIFY) ---
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

// Benchmark parameters
int batch_size;
int input_features;
int output_features;

// Input data
float* x;           // Input features: batch_size x input_features
float* h_prev;      // Previous hidden state: batch_size x output_features

// GRU weights and biases
float* Wz; float* Wr; float* Wh; // Input-to-hidden weights: input_features x output_features
float* Uz; float* Ur; float* Uh; // Hidden-to-hidden weights: output_features x output_features
float* bz; float* br; float* bh; // Biases: output_features

// Intermediate and output buffers
float* z_gate_input;
float* r_gate_input;
float* h_tilde_input;
float* z_gate;
float* r_gate;
float* h_tilde;
float* r_h_prod;
float* h_next;

// Final result accumulator
float final_result = 0.0f;

// Helper to generate a random float between -1.0 and 1.0
static float random_float() {
    return (mt_rand() / (float)UINT32_MAX) * 2.0f - 1.0f;
}

// Sigmoid activation function
static inline float sigmoid(float val) {
    return 1.0f / (1.0f + expf(-val));
}

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <batch_size> <input_features> <output_features> <seed>\n", argv[0]);
        exit(1);
    }

    batch_size = atoi(argv[1]);
    input_features = atoi(argv[2]);
    output_features = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    mt_seed(seed);

    // Allocate memory for inputs
    x = (float*)malloc(batch_size * input_features * sizeof(float));
    h_prev = (float*)malloc(batch_size * output_features * sizeof(float));

    // Allocate memory for weights and biases
    Wz = (float*)malloc(input_features * output_features * sizeof(float));
    Wr = (float*)malloc(input_features * output_features * sizeof(float));
    Wh = (float*)malloc(input_features * output_features * sizeof(float));
    Uz = (float*)malloc(output_features * output_features * sizeof(float));
    Ur = (float*)malloc(output_features * output_features * sizeof(float));
    Uh = (float*)malloc(output_features * output_features * sizeof(float));
    bz = (float*)malloc(output_features * sizeof(float));
    br = (float*)malloc(output_features * sizeof(float));
    bh = (float*)malloc(output_features * sizeof(float));
    
    // Allocate memory for intermediate and output buffers
    size_t out_buf_size = batch_size * output_features * sizeof(float);
    z_gate_input = (float*)malloc(out_buf_size);
    r_gate_input = (float*)malloc(out_buf_size);
    h_tilde_input = (float*)malloc(out_buf_size);
    z_gate = (float*)malloc(out_buf_size);
    r_gate = (float*)malloc(out_buf_size);
    h_tilde = (float*)malloc(out_buf_size);
    r_h_prod = (float*)malloc(out_buf_size);
    h_next = (float*)malloc(out_buf_size);

    // Initialize data with random values
    for(int i = 0; i < batch_size * input_features; ++i) x[i] = random_float();
    for(int i = 0; i < batch_size * output_features; ++i) h_prev[i] = random_float();
    for(int i = 0; i < input_features * output_features; ++i) Wz[i] = random_float();
    for(int i = 0; i < input_features * output_features; ++i) Wr[i] = random_float();
    for(int i = 0; i < input_features * output_features; ++i) Wh[i] = random_float();
    for(int i = 0; i < output_features * output_features; ++i) Uz[i] = random_float();
    for(int i = 0; i < output_features * output_features; ++i) Ur[i] = random_float();
    for(int i = 0; i < output_features * output_features; ++i) Uh[i] = random_float();
    for(int i = 0; i < output_features; ++i) bz[i] = random_float();
    for(int i = 0; i < output_features; ++i) br[i] = random_float();
    for(int i = 0; i < output_features; ++i) bh[i] = random_float();
}

void run_computation() {
    // --- Update Gate (z) --- 
    // z_gate_input = x @ Wz + h_prev @ Uz + bz

    // x @ Wz
    for (int i = 0; i < batch_size; ++i) {
        for (int j = 0; j < output_features; ++j) {
            float sum = 0.0f;
            for (int k = 0; k < input_features; ++k) {
                sum += x[i * input_features + k] * Wz[k * output_features + j];
            }
            z_gate_input[i * output_features + j] = sum;
        }
    }
    // += h_prev @ Uz
    for (int i = 0; i < batch_size; ++i) {
        for (int j = 0; j < output_features; ++j) {
            float sum = 0.0f;
            for (int k = 0; k < output_features; ++k) {
                sum += h_prev[i * output_features + k] * Uz[k * output_features + j];
            }
            z_gate_input[i * output_features + j] += sum;
        }
    }
    // += bz (broadcast) and apply sigmoid
    for (int i = 0; i < batch_size; ++i) {
        for (int j = 0; j < output_features; ++j) {
            z_gate[i * output_features + j] = sigmoid(z_gate_input[i * output_features + j] + bz[j]);
        }
    }

    // --- Reset Gate (r) ---
    // r_gate_input = x @ Wr + h_prev @ Ur + br
    
    // x @ Wr
    for (int i = 0; i < batch_size; ++i) {
        for (int j = 0; j < output_features; ++j) {
            float sum = 0.0f;
            for (int k = 0; k < input_features; ++k) {
                sum += x[i * input_features + k] * Wr[k * output_features + j];
            }
            r_gate_input[i * output_features + j] = sum;
        }
    }
    // += h_prev @ Ur
    for (int i = 0; i < batch_size; ++i) {
        for (int j = 0; j < output_features; ++j) {
            float sum = 0.0f;
            for (int k = 0; k < output_features; ++k) {
                sum += h_prev[i * output_features + k] * Ur[k * output_features + j];
            }
            r_gate_input[i * output_features + j] += sum;
        }
    }
    // += br (broadcast) and apply sigmoid
    for (int i = 0; i < batch_size; ++i) {
        for (int j = 0; j < output_features; ++j) {
            r_gate[i * output_features + j] = sigmoid(r_gate_input[i * output_features + j] + br[j]);
        }
    }

    // --- Candidate Hidden State (h_tilde) ---
    // h_tilde = tanh(x @ Wh + (r_gate .* h_prev) @ Uh + bh)

    // r_h_prod = r_gate .* h_prev (element-wise product)
    for (int i = 0; i < batch_size * output_features; ++i) {
        r_h_prod[i] = r_gate[i] * h_prev[i];
    }
    // x @ Wh
    for (int i = 0; i < batch_size; ++i) {
        for (int j = 0; j < output_features; ++j) {
            float sum = 0.0f;
            for (int k = 0; k < input_features; ++k) {
                sum += x[i * input_features + k] * Wh[k * output_features + j];
            }
            h_tilde_input[i * output_features + j] = sum;
        }
    }
    // += r_h_prod @ Uh
    for (int i = 0; i < batch_size; ++i) {
        for (int j = 0; j < output_features; ++j) {
            float sum = 0.0f;
            for (int k = 0; k < output_features; ++k) {
                sum += r_h_prod[i * output_features + k] * Uh[k * output_features + j];
            }
            h_tilde_input[i * output_features + j] += sum;
        }
    }
    // += bh (broadcast) and apply tanh
    for (int i = 0; i < batch_size; ++i) {
        for (int j = 0; j < output_features; ++j) {
            h_tilde[i * output_features + j] = tanhf(h_tilde_input[i * output_features + j] + bh[j]);
        }
    }

    // --- Final Hidden State (h_next) ---
    // h_next = (1 - z_gate) .* h_prev + z_gate .* h_tilde
    for (int i = 0; i < batch_size * output_features; ++i) {
        h_next[i] = (1.0f - z_gate[i]) * h_prev[i] + z_gate[i] * h_tilde[i];
    }

    // Accumulate the final result to prevent dead code elimination
    double acc = 0.0;
    for (int i = 0; i < batch_size * output_features; ++i) {
        acc += h_next[i];
    }
    final_result = (float) acc;
}

void cleanup() {
    free(x);
    free(h_prev);
    free(Wz); free(Wr); free(Wh);
    free(Uz); free(Ur); free(Uh);
    free(bz); free(br); free(bh);
    free(z_gate_input);
    free(r_gate_input);
    free(h_tilde_input);
    free(z_gate);
    free(r_gate);
    free(h_tilde);
    free(r_h_prod);
    free(h_next);
}

int main(int argc, char* argv[]) {
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