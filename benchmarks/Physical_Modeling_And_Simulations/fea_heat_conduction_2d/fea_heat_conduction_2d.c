#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- Mersenne Twister (MT19937) a_rand-in-c.h --- (Do Not Modify)
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

// Benchmark parameters and data structures
static int g_mesh_size_x;
static int g_mesh_size_y;
static int g_num_time_steps;
static int g_padded_x;
static int g_padded_y;
static double* g_temperature = NULL;
static double* g_temperature_next = NULL;
static double g_result = 0.0;

// The heat diffusion constant
const double ALPHA = 0.1;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s mesh_size_x mesh_size_y num_time_steps seed\n", argv[0]);
        exit(1);
    }

    g_mesh_size_x = atoi(argv[1]);
    g_mesh_size_y = atoi(argv[2]);
    g_num_time_steps = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    
    mt_seed(seed);

    if (g_mesh_size_x <= 0 || g_mesh_size_y <= 0 || g_num_time_steps <= 0) {
        fprintf(stderr, "FATAL: Mesh dimensions and time steps must be positive integers.\n");
        exit(1);
    }

    // Add padding for boundary conditions
    g_padded_x = g_mesh_size_x + 2;
    g_padded_y = g_mesh_size_y + 2;
    size_t grid_size_bytes = (size_t)g_padded_x * g_padded_y * sizeof(double);

    g_temperature = (double*)malloc(grid_size_bytes);
    g_temperature_next = (double*)malloc(grid_size_bytes);

    if (!g_temperature || !g_temperature_next) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize grids
    for (int y = 0; y < g_padded_y; ++y) {
        for (int x = 0; x < g_padded_x; ++x) {
            size_t idx = (size_t)y * g_padded_x + x;
            g_temperature[idx] = 0.0;
            g_temperature_next[idx] = 0.0;
        }
    }

    // Create a central hot-spot
    int hot_spot_x_start = g_mesh_size_x / 4;
    int hot_spot_x_end = g_mesh_size_x * 3 / 4;
    int hot_spot_y_start = g_mesh_size_y / 4;
    int hot_spot_y_end = g_mesh_size_y * 3 / 4;
    
    for (int y = hot_spot_y_start; y < hot_spot_y_end; ++y) {
        for (int x = hot_spot_x_start; x < hot_spot_x_end; ++x) {
            // Add 1 to account for padding
            size_t idx = (size_t)(y + 1) * g_padded_x + (x + 1);
            // Initialize with a high temperature plus minor random variations
            g_temperature[idx] = 1000.0 + (mt_rand() / (double)UINT32_MAX) * 10.0;
        }
    }
}

void run_computation() {
    double *current_temp = g_temperature;
    double *next_temp = g_temperature_next;

    for (int t = 0; t < g_num_time_steps; ++t) {
        // Apply 5-point stencil for heat diffusion
        for (int y = 1; y <= g_mesh_size_y; ++y) {
            for (int x = 1; x <= g_mesh_size_x; ++x) {
                size_t idx = (size_t)y * g_padded_x + x;
                double up = current_temp[idx - g_padded_x];
                double down = current_temp[idx + g_padded_x];
                double left = current_temp[idx - 1];
                double right = current_temp[idx + 1];
                double center = current_temp[idx];
                next_temp[idx] = center + ALPHA * (left + right + up + down - 4.0 * center);
            }
        }

        // Swap pointers for the next iteration
        double *swap = current_temp;
        current_temp = next_temp;
        next_temp = swap;
    }

    // To prevent DCE, calculate a checksum of the final temperatures
    double total_temp = 0.0;
    for (int y = 1; y <= g_mesh_size_y; ++y) {
        for (int x = 1; x <= g_mesh_size_x; ++x) {
            total_temp += current_temp[(size_t)y * g_padded_x + x];
        }
    }
    g_result = total_temp;
}

void cleanup() {
    free(g_temperature);
    free(g_temperature_next);
    g_temperature = NULL;
    g_temperature_next = NULL;
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
    printf("%f\n", g_result);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
