#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>
#include <string.h>

// --- MERSENNE TWISTER (DO NOT MODIFY) ---
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

// --- BENCHMARK DATA AND PARAMETERS ---
static int num_samples;
static int num_channels;
static int num_iterations;

// Input data: pre-processed (centered) mixed signals
// Dims: [num_channels][num_samples]
static float** mixed_signals;

// Output: the computed unmixing matrix
// Dims: [num_channels][num_channels]
static float** W;

// Final result to prevent dead code elimination
static double result_checksum;

// --- HELPER FUNCTIONS ---

// Allocate a 2D float matrix
float** allocate_matrix(int rows, int cols) {
    float** mat = (float**)malloc(rows * sizeof(float*));
    if (!mat) return NULL;
    for (int i = 0; i < rows; i++) {
        mat[i] = (float*)malloc(cols * sizeof(float));
        if (!mat[i]) {
            // Rollback allocation on failure
            for (int k = 0; k < i; k++) free(mat[k]);
            free(mat);
            return NULL;
        }
    }
    return mat;
}

// Free a 2D float matrix
void free_matrix(float** mat, int rows) {
    if (!mat) return;
    for (int i = 0; i < rows; i++) {
        free(mat[i]);
    }
    free(mat);
}

// Generate a random float in [-1.0, 1.0]
float rand_float() {
    return ((float)mt_rand() / (float)UINT32_MAX) * 2.0f - 1.0f;
}

// Non-quadratic non-linearity function g(u)
static inline float g(float u) {
    return tanhf(u);
}

// Derivative of g(u)
static inline float g_prime(float u) {
    float tanh_u = tanhf(u);
    return 1.0f - tanh_u * tanh_u;
}

// --- BENCHMARK FUNCTIONS --- 

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_samples num_channels num_iterations seed\n", argv[0]);
        exit(1);
    }

    num_samples = atoi(argv[1]);
    num_channels = atoi(argv[2]);
    num_iterations = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    if (num_samples <= 0 || num_channels <= 0 || num_iterations <= 0) {
        fprintf(stderr, "FATAL: a parameter is non-positive.\n");
        exit(1);
    }

    mt_seed(seed);

    mixed_signals = allocate_matrix(num_channels, num_samples);
    W = allocate_matrix(num_channels, num_channels);
    if (!mixed_signals || !W) {
        fprintf(stderr, "FATAL: Memory allocation failed\n");
        exit(1);
    }
    
    // Generate pre-processed (whitened) data
    // For this benchmark, we use centered random data as a substitute.
    for (int c = 0; c < num_channels; c++) {
        double sum = 0.0;
        for (int s = 0; s < num_samples; s++) {
            mixed_signals[c][s] = rand_float();
            sum += mixed_signals[c][s];
        }
        float mean = (float)(sum / num_samples);
        for (int s = 0; s < num_samples; s++) {
            mixed_signals[c][s] -= mean;
        }
    }
}

void run_computation() {
    float* w = (float*)malloc(num_channels * sizeof(float));
    float* w_new = (float*)malloc(num_channels * sizeof(float));

    // Deflationary approach: find one component at a time
    for (int p = 0; p < num_channels; p++) {
        
        // Initialize weight vector w with random values and normalize
        float norm = 0.0f;
        for (int i = 0; i < num_channels; i++) {
            w[i] = rand_float();
            norm += w[i] * w[i];
        }
        norm = sqrtf(norm);
        for (int i = 0; i < num_channels; i++) {
            w[i] /= norm;
        }

        // Main FastICA iteration loop for one component
        for (int iter = 0; iter < num_iterations; iter++) {
            // 1. Compute expectations for the update rule
            memset(w_new, 0, num_channels * sizeof(float));
            float mean_g_prime_wTx = 0.0f;

            for (int s = 0; s < num_samples; s++) {
                float wTx = 0.0f;
                // Dot product: w_transpose * x_s
                for (int c = 0; c < num_channels; c++) {
                    wTx += w[c] * mixed_signals[c][s];
                }
                float g_wTx = g(wTx);
                mean_g_prime_wTx += g_prime(wTx);
                // Accumulate E{x * g(w'x)}
                for (int c = 0; c < num_channels; c++) {
                    w_new[c] += mixed_signals[c][s] * g_wTx;
                }
            }
            mean_g_prime_wTx /= num_samples;
            for (int c = 0; c < num_channels; c++) {
                w_new[c] /= num_samples;
            }
            
            // 2. Apply the FastICA update rule
            // w_new = E{x * g(w'x)} - E{g'(w'x)} * w
            for (int c = 0; c < num_channels; c++) {
                w_new[c] -= mean_g_prime_wTx * w[c];
            }

            // 3. Decorrelate w_new from previously found components (Gram-Schmidt)
            for (int j = 0; j < p; j++) {
                float proj = 0.0f;
                for (int c = 0; c < num_channels; c++) {
                    proj += w_new[c] * W[j][c];
                }
                for (int c = 0; c < num_channels; c++) {
                    w_new[c] -= proj * W[j][c];
                }
            }

            // 4. Normalize the new weight vector
            float norm_new = 0.0f;
            for (int c = 0; c < num_channels; c++) {
                norm_new += w_new[c] * w_new[c];
            }
            norm_new = sqrtf(norm_new);
            if (norm_new > 1e-9f) { // Avoid division by zero
                for (int c = 0; c < num_channels; c++) {
                    w[c] = w_new[c] / norm_new;
                }
            }
        }
        // Store the found component vector in the unmixing matrix W
        memcpy(W[p], w, num_channels * sizeof(float));
    }

    // Calculate a checksum of the result matrix to prevent dead-code elimination
    result_checksum = 0.0;
    for (int i = 0; i < num_channels; i++) {
        for (int j = 0; j < num_channels; j++) {
            result_checksum += W[i][j];
        }
    }

    free(w);
    free(w_new);
}

void cleanup() {
    free_matrix(mixed_signals, num_channels);
    free_matrix(W, num_channels);
}

// --- MAIN FUNCTION --- 

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%.6f\n", result_checksum);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
