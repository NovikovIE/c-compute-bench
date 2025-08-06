#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <limits.h>

// MERSENNE TWISTER (DO NOT MODIFY)
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

// --- BENCHMARK DATA AND PARAMETERS ---
int grid_width;
int grid_height;
int kernel_size;

float* input_grid;
float* output_grid;
double final_result;

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s grid_width grid_height kernel_size seed\n", argv[0]);
        exit(1);
    }

    grid_width = atoi(argv[1]);
    grid_height = atoi(argv[2]);
    kernel_size = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    if (grid_width <= 0 || grid_height <= 0 || kernel_size <= 0) {
        fprintf(stderr, "Error: grid dimensions and kernel size must be positive.\n");
        exit(1);
    }
    if (kernel_size % 2 == 0) {
        fprintf(stderr, "Error: kernel_size must be an odd number.\n");
        exit(1);
    }
    if (kernel_size > grid_width || kernel_size > grid_height) {
         fprintf(stderr, "Error: kernel is larger than the grid dimensions.\n");
         exit(1);
    }

    mt_seed(seed);

    size_t num_elements = (size_t)grid_width * grid_height;
    input_grid = (float*)malloc(num_elements * sizeof(float));
    output_grid = (float*)calloc(num_elements, sizeof(float)); // Initializes to 0.0

    if (!input_grid || !output_grid) {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        exit(1);
    }

    for (size_t i = 0; i < num_elements; i++) {
        input_grid[i] = (float)mt_rand() / (float)UINT32_MAX;
    }
}

void run_computation() {
    int kernel_radius = kernel_size / 2;
    double kernel_divisor = (double)kernel_size * kernel_size;

    // Iterate over the grid, avoiding the borders where the kernel would go out of bounds
    for (int y = kernel_radius; y < grid_height - kernel_radius; ++y) {
        for (int x = kernel_radius; x < grid_width - kernel_radius; ++x) {
            double sum = 0.0;

            // Apply the kernel by summing up neighbors
            for (int ky = -kernel_radius; ky <= kernel_radius; ++ky) {
                for (int kx = -kernel_radius; kx <= kernel_radius; ++kx) {
                    int sample_y = y + ky;
                    int sample_x = x + kx;
                    sum += input_grid[(size_t)sample_y * grid_width + sample_x];
                }
            }
            output_grid[(size_t)y * grid_width + x] = (float)(sum / kernel_divisor);
        }
    }

    // Calculate a final result to prevent dead code elimination.
    // Summing all elements is fine as uncomputed border cells are 0.0 from calloc.
    final_result = 0.0;
    size_t num_elements = (size_t)grid_width * grid_height;
    for (size_t i = 0; i < num_elements; i++) {
        final_result += output_grid[i];
    }
}

void cleanup() {
    free(input_grid);
    free(output_grid);
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

    // Print final result to stdout
    printf("%f\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
