#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- BEGIN MERSENNE TWISTER --- (DO NOT MODIFY)
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

// Benchmark parameters
int grid_size_x;
int grid_size_y;
int time_steps;

// Data grids
double** grid;
double** next_grid;

// Result checksum
double result_checksum;

// Function to allocate a 2D grid
double** allocate_grid(int nx, int ny) {
    double** g = (double**)malloc(nx * sizeof(double*));
    if (g == NULL) return NULL;

    for (int i = 0; i < nx; ++i) {
        g[i] = (double*)malloc(ny * sizeof(double));
        if (g[i] == NULL) {
            // Cleanup already allocated rows
            for (int k = 0; k < i; ++k) free(g[k]);
            free(g);
            return NULL;
        }
    }
    return g;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s grid_size_x grid_size_y time_steps seed\n", argv[0]);
        exit(1);
    }

    grid_size_x = atoi(argv[1]);
    grid_size_y = atoi(argv[2]);
    time_steps = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    
    mt_seed(seed);

    grid = allocate_grid(grid_size_x, grid_size_y);
    next_grid = allocate_grid(grid_size_x, grid_size_y);

    if (grid == NULL || next_grid == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize grids
    // Set initial conditions randomly, with boundaries held at 0.
    for (int i = 0; i < grid_size_x; ++i) {
        for (int j = 0; j < grid_size_y; ++j) {
            if (i == 0 || i == grid_size_x - 1 || j == 0 || j == grid_size_y - 1) {
                grid[i][j] = 0.0;
                next_grid[i][j] = 0.0;
            } else {
                grid[i][j] = (double)mt_rand() / (double)UINT32_MAX; // Random value between 0.0 and 1.0
                next_grid[i][j] = 0.0;
            }
        }
    }
    result_checksum = 0.0;
}

void run_computation() {
    // Simulate the 2D heat equation using a 5-point stencil finite difference method.
    // PDE: dU/dt = alpha * (d^2U/dx^2 + d^2U/dy^2)
    const double alpha = 0.1; // Diffusion coefficient

    for (int t = 0; t < time_steps; ++t) {
        // Update interior points
        for (int i = 1; i < grid_size_x - 1; ++i) {
            for (int j = 1; j < grid_size_y - 1; ++j) {
                next_grid[i][j] = grid[i][j] + alpha * (
                    grid[i+1][j] + 
                    grid[i-1][j] + 
                    grid[i][j+1] + 
                    grid[i][j-1] - 
                    4.0 * grid[i][j]
                );
            }
        }

        // Swap grids for the next iteration
        double** temp = grid;
        grid = next_grid;
        next_grid = temp;
    }

    // Calculate a checksum to prevent dead code elimination
    for (int i = 0; i < grid_size_x; ++i) {
        for (int j = 0; j < grid_size_y; ++j) {
            result_checksum += grid[i][j];
        }
    }
}

void cleanup() {
    for (int i = 0; i < grid_size_x; ++i) {
        free(grid[i]);
        free(next_grid[i]);
    }
    free(grid);
    free(next_grid);
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
    printf("%f\n", result_checksum);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
