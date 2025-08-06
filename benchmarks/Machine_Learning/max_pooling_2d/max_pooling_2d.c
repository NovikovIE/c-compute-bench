/**
 * @file max_pooling_2d.c
 * @brief A benchmark for the 2D Max Pooling operation in machine learning.
 *
 * This program implements and times the 2D max pooling operation, a common
 * downsampling strategy used in Convolutional Neural Networks (CNNs). It operates
 * on a 3D input tensor (height, width, channels) and produces a smaller 3D output
 * tensor by taking the maximum value within non-overlapping pooling windows.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h> // For -INFINITY

// --- Mersenne Twister (MT19937) --- (DO NOT MODIFY)
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
// --- End of Mersenne Twister ---

// Benchmark parameters
static int INPUT_HEIGHT, INPUT_WIDTH, CHANNELS, POOL_SIZE;
static int OUTPUT_HEIGHT, OUTPUT_WIDTH;

// Data buffers
static float *input_tensor;
static float *output_tensor;

// Final result accumulator
static double final_result = 0.0;

/**
 * @brief Parses arguments, allocates memory, and initializes input data.
 * 
 * Allocates the input and output tensors on the heap and fills the input
 * tensor with random floating-point values between 0.0 and 1.0.
 */
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <input_height> <input_width> <channels> <pool_size> <seed>\n", argv[0]);
        exit(1);
    }

    INPUT_HEIGHT = atoi(argv[1]);
    INPUT_WIDTH = atoi(argv[2]);
    CHANNELS = atoi(argv[3]);
    POOL_SIZE = atoi(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    if (INPUT_HEIGHT <= 0 || INPUT_WIDTH <= 0 || CHANNELS <= 0 || POOL_SIZE <= 0) {
        fprintf(stderr, "ERROR: All parameters must be positive integers.\n");
        exit(1);
    }
    if (INPUT_HEIGHT % POOL_SIZE != 0 || INPUT_WIDTH % POOL_SIZE != 0) {
        fprintf(stderr, "ERROR: Input dimensions must be divisible by pool_size.\n");
        exit(1);
    }

    mt_seed(seed);

    OUTPUT_HEIGHT = INPUT_HEIGHT / POOL_SIZE;
    OUTPUT_WIDTH = INPUT_WIDTH / POOL_SIZE;

    size_t input_size = (size_t)INPUT_HEIGHT * INPUT_WIDTH * CHANNELS;
    size_t output_size = (size_t)OUTPUT_HEIGHT * OUTPUT_WIDTH * CHANNELS;

    input_tensor = (float *)malloc(input_size * sizeof(float));
    output_tensor = (float *)malloc(output_size * sizeof(float));

    if (!input_tensor || !output_tensor) {
        fprintf(stderr, "ERROR: Memory allocation failed.\n");
        exit(1);
    }

    for (size_t i = 0; i < input_size; i++) {
        input_tensor[i] = (float)mt_rand() / (float)UINT32_MAX;
    }
}

/**
 * @brief Executes the core max pooling computation.
 *
 * Iterates through the input tensor, applying a max pooling filter to each
 * window and storing the result in the output tensor. Finally, it computes a
 * sum of all output values to serve as a verifiable result and prevent dead
 * code elimination by the compiler.
 */
void run_computation() {
    long input_slice_size = (long)INPUT_HEIGHT * INPUT_WIDTH;
    long output_slice_size = (long)OUTPUT_HEIGHT * OUTPUT_WIDTH;

    for (int c = 0; c < CHANNELS; c++) {
        for (int oh = 0; oh < OUTPUT_HEIGHT; oh++) {
            for (int ow = 0; ow < OUTPUT_WIDTH; ow++) {
                float max_val = -INFINITY;
                int start_h = oh * POOL_SIZE;
                int start_w = ow * POOL_SIZE;

                for (int ph = 0; ph < POOL_SIZE; ph++) {
                    for (int pw = 0; pw < POOL_SIZE; pw++) {
                        int ih = start_h + ph;
                        int iw = start_w + pw;
                        long input_idx = c * input_slice_size + (long)ih * INPUT_WIDTH + iw;
                        if (input_tensor[input_idx] > max_val) {
                            max_val = input_tensor[input_idx];
                        }
                    }
                }
                long output_idx = c * output_slice_size + (long)oh * OUTPUT_WIDTH + ow;
                output_tensor[output_idx] = max_val;
            }
        }
    }

    // Accumulate the result to prevent dead code elimination
    double sum = 0.0;
    size_t output_size = (size_t)OUTPUT_HEIGHT * OUTPUT_WIDTH * CHANNELS;
    for (size_t i = 0; i < output_size; i++) {
        sum += output_tensor[i];
    }
    final_result = sum;
}

/**
 * @brief Frees all memory allocated in setup_benchmark.
 */
void cleanup() {
    free(input_tensor);
    free(output_tensor);
}

/**
 * @brief Main function to orchestrate the benchmark.
 *
 * This function follows a strict sequence:
 * 1. Setup: Prepare data for the benchmark.
 * 2. Time Computation: Measure the execution time of the core algorithm.
 * 3. Cleanup: Release all allocated resources.
 * 4. Output: Print the final result to stdout and the time to stderr.
 */
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final result to stdout
    printf("%f\n", final_result);

    // Print the execution time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
