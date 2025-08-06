#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- Mersenne Twister (MT19937) Generator - Do Not Modify ---
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

// --- Benchmark Globals ---
int image_height, image_width, input_channels, output_channels, kernel_size, stride;
int output_height, output_width;

float *input_image; 
float *kernels;     
float *output_image;

float final_result = 0.0f;

// --- Benchmark Functions ---

// Function to generate a random float between -1.0 and 1.0
float rand_float() {
    return ((float)mt_rand() / (float)(UINT32_MAX / 2)) - 1.0f;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 8) {
        fprintf(stderr, "Usage: %s image_height image_width input_channels output_channels kernel_size stride seed\n", argv[0]);
        exit(1);
    }

    image_height = atoi(argv[1]);
    image_width = atoi(argv[2]);
    input_channels = atoi(argv[3]);
    output_channels = atoi(argv[4]);
    kernel_size = atoi(argv[5]);
    stride = atoi(argv[6]);
    uint32_t seed = (uint32_t)atoi(argv[7]);

    mt_seed(seed);

    // Calculate output dimensions
    output_height = (image_height - kernel_size) / stride + 1;
    output_width = (image_width - kernel_size) / stride + 1;

    // Allocate memory
    unsigned long long input_size = (unsigned long long)input_channels * image_height * image_width;
    unsigned long long kernel_size_full = (unsigned long long)output_channels * input_channels * kernel_size * kernel_size;
    unsigned long long output_size = (unsigned long long)output_channels * output_height * output_width;

    input_image = (float*)malloc(input_size * sizeof(float));
    kernels = (float*)malloc(kernel_size_full * sizeof(float));
    output_image = (float*)malloc(output_size * sizeof(float));

    if (!input_image || !kernels || !output_image) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize input image and kernels with random data
    for (unsigned long long i = 0; i < input_size; ++i) {
        input_image[i] = rand_float();
    }

    for (unsigned long long i = 0; i < kernel_size_full; ++i) {
        kernels[i] = rand_float();
    }
}

void run_computation() {
    // For each output feature map
    for (int oc = 0; oc < output_channels; ++oc) {
        // For each pixel in the output feature map
        for (int oy = 0; oy < output_height; ++oy) {
            for (int ox = 0; ox < output_width; ++ox) {
                float sum = 0.0f;
                int iy_base = oy * stride;
                int ix_base = ox * stride;

                // For each input channel
                for (int ic = 0; ic < input_channels; ++ic) {
                    // Apply the kernel
                    for (int ky = 0; ky < kernel_size; ++ky) {
                        for (int kx = 0; kx < kernel_size; ++kx) {
                            int iy = iy_base + ky;
                            int ix = ix_base + kx;

                            // Index calculations for flattened arrays
                            long long input_idx = (long long)ic * image_height * image_width + (long long)iy * image_width + ix;
                            long long kernel_idx = (long long)oc * input_channels * kernel_size * kernel_size + 
                                                 (long long)ic * kernel_size * kernel_size + 
                                                 (long long)ky * kernel_size + kx;

                            sum += input_image[input_idx] * kernels[kernel_idx];
                        }
                    }
                }
                output_image[(long long)oc * output_height * output_width + (long long)oy * output_width + ox] = sum;
            }
        }
    }

    // Accumulate a result to prevent dead code elimination
    float total_sum = 0.0f;
    unsigned long long output_size = (unsigned long long)output_channels * output_height * output_width;
    for (unsigned long long i = 0; i < output_size; ++i) {
        total_sum += output_image[i];
    }
    final_result = total_sum;
}

void cleanup() {
    free(input_image);
    free(kernels);
    free(output_image);
    input_image = NULL;
    kernels = NULL;
    output_image = NULL;
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
