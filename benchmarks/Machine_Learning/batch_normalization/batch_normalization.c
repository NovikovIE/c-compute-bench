#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- START of Mersenne Twister (Do Not Modify) ---
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
// --- END of Mersenne Twister ---

// --- Benchmark Globals ---
int BATCH_SIZE;
int NUM_FEATURES;
const float EPSILON = 1e-5f;

// Data pointers for contiguous 2D arrays
float** input_batch;
float*  input_batch_data;
float** output_batch;
float*  output_batch_data;

// Learnable parameters and buffers
float* bn_gamma;
float* beta;
float* mean;
float* variance;

// Result accumulator
double final_result;

// Helper function to generate a random float between 0.0 and 1.0
float rand_float() {
    return (float)mt_rand() / (float)UINT32_MAX;
}

// Helper to allocate a contiguous 2D float array and its pointer array
float** allocate_matrix(int rows, int cols, float** data_blob) {
    float** matrix = (float**)malloc(rows * sizeof(float*));
    if (!matrix) return NULL;
    
    *data_blob = (float*)malloc((size_t)rows * cols * sizeof(float));
    if (!*data_blob) {
        free(matrix);
        return NULL;
    }

    for (int i = 0; i < rows; ++i) {
        matrix[i] = &((*data_blob)[(size_t)i * cols]);
    }
    return matrix;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_features> <batch_size> <seed>\n", argv[0]);
        exit(1);
    }

    NUM_FEATURES = atoi(argv[1]);
    BATCH_SIZE = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);
    mt_seed(seed);

    if (NUM_FEATURES <= 0 || BATCH_SIZE <= 0) {
        fprintf(stderr, "FATAL: num_features and batch_size must be positive integers.\n");
        exit(1);
    }

    // Allocate memories
    input_batch = allocate_matrix(BATCH_SIZE, NUM_FEATURES, &input_batch_data);
    output_batch = allocate_matrix(BATCH_SIZE, NUM_FEATURES, &output_batch_data);
    bn_gamma = (float*)malloc(NUM_FEATURES * sizeof(float));
    beta = (float*)malloc(NUM_FEATURES * sizeof(float));
    mean = (float*)malloc(NUM_FEATURES * sizeof(float));
    variance = (float*)malloc(NUM_FEATURES * sizeof(float));

    if (!input_batch || !output_batch || !bn_gamma || !beta || !mean || !variance) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize input data and parameters
    for (int i = 0; i < BATCH_SIZE; ++i) {
        for (int j = 0; j < NUM_FEATURES; ++j) {
            input_batch[i][j] = rand_float() * 10.0f - 5.0f; // Random data in [-5, 5]
        }
    }

    for (int j = 0; j < NUM_FEATURES; ++j) {
        bn_gamma[j] = rand_float() + 0.5f;   // Scaling parameter, away from zero
        beta[j] = rand_float() - 0.5f;    // Shifting parameter
    }
    
    final_result = 0.0;
}

void run_computation() {
    // --- Batch Normalization Forward Pass ---

    double batch_size_inv = 1.0 / BATCH_SIZE;

    // Step 1 & 2: Calculate mean and variance for each feature
    for (int j = 0; j < NUM_FEATURES; ++j) {
        double sum = 0.0;
        for (int i = 0; i < BATCH_SIZE; ++i) {
            sum += input_batch[i][j];
        }
        mean[j] = (float)(sum * batch_size_inv);

        double var_sum = 0.0;
        for (int i = 0; i < BATCH_SIZE; ++i) {
            float diff = input_batch[i][j] - mean[j];
            var_sum += diff * diff;
        }
        variance[j] = (float)(var_sum * batch_size_inv);
    }

    // Step 3: Normalize, scale, and shift
    for (int j = 0; j < NUM_FEATURES; ++j) {
        float inv_stddev = 1.0f / sqrtf(variance[j] + EPSILON);
        for (int i = 0; i < BATCH_SIZE; ++i) {
            float normalized = (input_batch[i][j] - mean[j]) * inv_stddev;
            output_batch[i][j] = normalized * bn_gamma[j] + beta[j];
        }
    }

    // --- Accumulate result to prevent dead code elimination ---
    double total_sum = 0.0;
    for (int i = 0; i < BATCH_SIZE; ++i) {
        for (int j = 0; j < NUM_FEATURES; ++j) {
            total_sum += output_batch[i][j];
        }
    }
    final_result = total_sum;
}

void cleanup() {
    free(input_batch_data);
    free(input_batch);
    free(output_batch_data);
    free(output_batch);
    free(bn_gamma);
    free(beta);
    free(mean);
    free(variance);
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
