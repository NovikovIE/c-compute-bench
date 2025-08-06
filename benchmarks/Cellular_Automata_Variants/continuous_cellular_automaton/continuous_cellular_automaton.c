#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- Mersenne Twister (MT19937) Generator - DO NOT MODIFY ---
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

// --- Benchmark Globals ---
int grid_size;
int num_steps;
int kernel_radius;

float *grid_a;
float *grid_b;
float *kernel;

float final_result;

// --- Benchmark Functions ---

// Generates a random float between 0.0 and 1.0
float random_float() {
    return (float)mt_rand() / (float)UINT32_MAX;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s grid_size num_steps kernel_radius seed\n", argv[0]);
        exit(1);
    }

    grid_size = atoi(argv[1]);
    num_steps = atoi(argv[2]);
    kernel_radius = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    mt_seed(seed);

    size_t grid_elements = (size_t)grid_size * grid_size;
    grid_a = (float *)malloc(grid_elements * sizeof(float));
    grid_b = (float *)malloc(grid_elements * sizeof(float));

    int kernel_dim = 2 * kernel_radius + 1;
    size_t kernel_elements = (size_t)kernel_dim * kernel_dim;
    kernel = (float *)malloc(kernel_elements * sizeof(float));

    if (!grid_a || !grid_b || !kernel) {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize grid with random float values
    for (size_t i = 0; i < grid_elements; ++i) {
        grid_a[i] = random_float();
    }

    // Initialize kernel with random weights (centered around 0)
    for (size_t i = 0; i < kernel_elements; ++i) {
        kernel[i] = random_float() - 0.5f;
    }
}

void run_computation() {
    float *read_grid = grid_a;
    float *write_grid = grid_b;

    int kernel_dim = 2 * kernel_radius + 1;
    float kernel_norm = 1.0f / (float)(kernel_dim * kernel_dim);

    for (int step = 0; step < num_steps; ++step) {
        for (int y = 0; y < grid_size; ++y) {
            for (int x = 0; x < grid_size; ++x) {
                float sum = 0.0f;

                // Apply convolution kernel
                for (int ky = -kernel_radius; ky <= kernel_radius; ++ky) {
                    for (int kx = -kernel_radius; kx <= kernel_radius; ++kx) {
                        // Toroidal (wrapping) boundary conditions
                        int neighbor_y = (y + ky + grid_size) % grid_size;
                        int neighbor_x = (x + kx + grid_size) % grid_size;

                        float neighbor_val = read_grid[neighbor_y * grid_size + neighbor_x];
                        float kernel_val = kernel[(ky + kernel_radius) * kernel_dim + (kx + kernel_radius)];
                        sum += neighbor_val * kernel_val;
                    }
                }

                // Apply a non-linear activation function (sigmoid-like)
                // The sum is normalized by kernel size to prevent extreme values
                float new_value = 1.0f / (1.0f + expf(-sum * kernel_norm));
                write_grid[y * grid_size + x] = new_value;
            }
        }

        // Swap grids for the next iteration
        float *temp = read_grid;
        read_grid = write_grid;
        write_grid = temp;
    }

    // Accumulate a final result to prevent dead code elimination
    double total_sum = 0.0;
    size_t grid_elements = (size_t)grid_size * grid_size;
    for (size_t i = 0; i < grid_elements; ++i) {
        total_sum += read_grid[i];
    }
    final_result = (float)total_sum;
}

void cleanup() {
    free(grid_a);
    free(grid_b);
    free(kernel);
}

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

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
