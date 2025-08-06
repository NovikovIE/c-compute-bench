#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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
// --- END MERSENNE TWISTER ---

// --- BENCHMARK DATA & GLOBALS ---
typedef struct {
    int num_dimensions;
    int grid_points_per_dim;
    double search_min;       // Lower bound of search space for each dimension
    double search_max;       // Upper bound of search space for each dimension
    double step_size;        // Distance between grid points
    double min_rastrigin_val; // Accumulated result: the minimum value found
} BenchmarkData;

static BenchmarkData g_data;

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_dimensions> <grid_points_per_dim> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_dimensions = atoi(argv[1]);
    g_data.grid_points_per_dim = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (g_data.num_dimensions <= 0 || g_data.grid_points_per_dim <= 1) {
        fprintf(stderr, "Error: num_dimensions must be > 0 and grid_points_per_dim > 1.\n");
        exit(1);
    }

    mt_seed(seed); // Seed generator as required, even if unused in this algorithm.

    // Rastrigin function is typically evaluated on the hypercube xi ∈ [-5.12, 5.12]
    g_data.search_min = -5.12;
    g_data.search_max = 5.12;
    g_data.step_size = (g_data.search_max - g_data.search_min) / (double)(g_data.grid_points_per_dim - 1);

    // Initialize the minimum value to positive infinity.
    g_data.min_rastrigin_val = 1.0/0.0;
}

// The Rastrigin function, which we aim to minimize.
// f(x) = 10n + ∑(xi^2 - 10cos(2πxi))
// Global minimum is f(0,0,...,0) = 0
double evaluate_rastrigin(double* point, int dimensions) {
    const double A = 10.0;
    double sum = A * dimensions;
    for (int i = 0; i < dimensions; i++) {
        double xi = point[i];
        sum += (xi * xi) - (A * cos(2.0 * M_PI * xi));
    }
    return sum;
}

void run_computation() {
    // Allocate memory for tracking grid position and coordinates.
    int* indices = (int*)calloc(g_data.num_dimensions, sizeof(int));
    if (!indices) { perror("Failed to allocate indices"); exit(1); }

    double* current_point = (double*)malloc(g_data.num_dimensions * sizeof(double));
    if (!current_point) { perror("Failed to allocate current_point"); free(indices); exit(1); }
    
    while (1) {
        // 1. Construct the current grid point's coordinates from indices.
        for (int i = 0; i < g_data.num_dimensions; i++) {
            current_point[i] = g_data.search_min + indices[i] * g_data.step_size;
        }

        // 2. Evaluate the function at this point.
        double current_val = evaluate_rastrigin(current_point, g_data.num_dimensions);

        // 3. Update the global minimum if the current value is smaller.
        if (current_val < g_data.min_rastrigin_val) {
            g_data.min_rastrigin_val = current_val;
        }

        // 4. Increment the multi-dimensional grid counter (like an odometer).
        int dim_to_inc = 0;
        while (dim_to_inc < g_data.num_dimensions) {
            indices[dim_to_inc]++;
            if (indices[dim_to_inc] < g_data.grid_points_per_dim) {
                break; // No carry-over, continue to the next grid point.
            }
            indices[dim_to_inc] = 0; // Reset current dimension and carry to the next.
            dim_to_inc++;
        }

        // If the carry-over propagates past the last dimension, we are done.
        if (dim_to_inc == g_data.num_dimensions) {
            break;
        }
    }

    free(indices);
    free(current_point);
}

void cleanup() {
    // No persistent heap allocations performed in setup_benchmark.
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

    // Print the final result (minimum value found) to stdout.
    printf("%f\n", g_data.min_rastrigin_val);

    // Print the execution time to stderr.
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
