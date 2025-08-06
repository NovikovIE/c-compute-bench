#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <string.h>

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

// Benchmark parameters and data structures
static int grid_size;
static int num_iterations;
static unsigned int final_result;

// Two grids for storing current and next states
static uint8_t** grid_a;
static uint8_t** grid_b;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <grid_size> <num_iterations> <seed>\n", argv[0]);
        exit(1);
    }

    grid_size = atoi(argv[1]);
    num_iterations = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (grid_size <= 0 || num_iterations <= 0) {
        fprintf(stderr, "FATAL: grid_size and num_iterations must be positive integers.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate memory for the grids (as 2D arrays)
    grid_a = (uint8_t**)malloc(grid_size * sizeof(uint8_t*));
    grid_b = (uint8_t**)malloc(grid_size * sizeof(uint8_t*));
    if (!grid_a || !grid_b) {
        fprintf(stderr, "FATAL: Memory allocation failed for grid pointers.\n");
        exit(1);
    }

    for (int i = 0; i < grid_size; ++i) {
        grid_a[i] = (uint8_t*)malloc(grid_size * sizeof(uint8_t));
        grid_b[i] = (uint8_t*)malloc(grid_size * sizeof(uint8_t));
        if (!grid_a[i] || !grid_b[i]) {
            fprintf(stderr, "FATAL: Memory allocation failed for grid rows.\n");
            exit(1);
        }
    }

    // Initialize grid_a with random states (0 or 1) and grid_b to zero
    for (int i = 0; i < grid_size; ++i) {
        for (int j = 0; j < grid_size; ++j) {
            grid_a[i][j] = mt_rand() % 2;
            grid_b[i][j] = 0;
        }
    }
}

void run_computation() {
    uint8_t** current_grid = grid_a;
    uint8_t** next_grid = grid_b;

    // Offsets for the 8 neighbors in a Moore neighborhood
    const int dx[] = {-1, -1, -1, 0, 0, 1, 1, 1};
    const int dy[] = {-1, 0, 1, -1, 1, -1, 0, 1};

    for (int iter = 0; iter < num_iterations; ++iter) {
        for (int i = 0; i < grid_size; ++i) {
            for (int j = 0; j < grid_size; ++j) {
                int live_neighbors = 0;

                // Count live neighbors with toroidal boundary conditions
                for (int k = 0; k < 8; ++k) {
                    int ni = (i + dx[k] + grid_size) % grid_size;
                    int nj = (j + dy[k] + grid_size) % grid_size;
                    live_neighbors += current_grid[ni][nj];
                }

                // Apply the Majority Vote rule
                if (live_neighbors > 4) {
                    next_grid[i][j] = 1;
                } else if (live_neighbors < 4) {
                    next_grid[i][j] = 0;
                } else { // if live_neighbors == 4, state is unchanged
                    next_grid[i][j] = current_grid[i][j];
                }
            }
        }

        // Swap grids for the next iteration
        uint8_t** temp = current_grid;
        current_grid = next_grid;
        next_grid = temp;
    }

    // Calculate a final result to prevent dead code elimination
    unsigned int sum = 0;
    for (int i = 0; i < grid_size; ++i) {
        for (int j = 0; j < grid_size; ++j) {
            sum += current_grid[i][j];
        }
    }
    final_result = sum;
}

void cleanup() {
    for (int i = 0; i < grid_size; ++i) {
        free(grid_a[i]);
        free(grid_b[i]);
    }
    free(grid_a);
    free(grid_b);
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
    printf("%u\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
