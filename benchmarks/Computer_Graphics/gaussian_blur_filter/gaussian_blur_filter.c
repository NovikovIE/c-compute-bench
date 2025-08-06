#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

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

// --- Benchmark Globals ---
int image_width;
int image_height;
int kernel_radius;

float *input_image = NULL;
float *output_image = NULL;
float *kernel = NULL;

double final_result = 0.0;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s image_width image_height kernel_radius seed\n", argv[0]);
        exit(1);
    }

    image_width = atoi(argv[1]);
    image_height = atoi(argv[2]);
    kernel_radius = atoi(argv[3]);
    uint32_t seed = (uint32_t)strtoul(argv[4], NULL, 10);

    mt_seed(seed);

    size_t image_size = (size_t)image_width * image_height;
    input_image = (float *)malloc(image_size * sizeof(float));
    output_image = (float *)malloc(image_size * sizeof(float));

    if (!input_image || !output_image) {
        fprintf(stderr, "Error: Memory allocation failed for images.\n");
        exit(1);
    }

    for (size_t i = 0; i < image_size; i++) {
        input_image[i] = (float)mt_rand() / (float)UINT32_MAX * 255.0f;
    }

    int kernel_dim = 2 * kernel_radius + 1;
    kernel = (float *)malloc(kernel_dim * kernel_dim * sizeof(float));
    if (!kernel) {
        fprintf(stderr, "Error: Memory allocation failed for kernel.\n");
        exit(1);
    }

    double sigma = (double)kernel_radius / 3.0;
    double s = 2.0 * sigma * sigma;
    double sum = 0.0;

    for (int y = -kernel_radius; y <= kernel_radius; y++) {
        for (int x = -kernel_radius; x <= kernel_radius; x++) {
            double r = sqrt(x * x + y * y);
            int kernel_idx = (y + kernel_radius) * kernel_dim + (x + kernel_radius);
            kernel[kernel_idx] = (float)(exp(-(r * r) / s) / (M_PI * s));
            sum += kernel[kernel_idx];
        }
    }

    // Normalize kernel
    for (int i = 0; i < kernel_dim * kernel_dim; i++) {
        kernel[i] /= sum;
    }
}

void run_computation() {
    int kernel_dim = 2 * kernel_radius + 1;
    for (int y = 0; y < image_height; y++) {
        for (int x = 0; x < image_width; x++) {
            float sum = 0.0f;
            for (int ky = -kernel_radius; ky <= kernel_radius; ky++) {
                for (int kx = -kernel_radius; kx <= kernel_radius; kx++) {
                    int src_x = x + kx;
                    int src_y = y + ky;

                    // Clamp to edge
                    if (src_x < 0) src_x = 0;
                    if (src_x >= image_width) src_x = image_width - 1;
                    if (src_y < 0) src_y = 0;
                    if (src_y >= image_height) src_y = image_height - 1;

                    float pixel_val = input_image[src_y * image_width + src_x];
                    float kernel_val = kernel[(ky + kernel_radius) * kernel_dim + (kx + kernel_radius)];
                    sum += pixel_val * kernel_val;
                }
            }
            output_image[y * image_width + x] = sum;
        }
    }

    // Calculate a checksum to prevent dead code elimination
    double checksum = 0.0;
    size_t image_size = (size_t)image_width * image_height;
    for (size_t i = 0; i < image_size; i++) {
        checksum += output_image[i];
    }
    final_result = checksum;
}

void cleanup() {
    if (input_image) free(input_image);
    if (output_image) free(output_image);
    if (kernel) free(kernel);
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
