#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <string.h>

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
static int grid_size;
static int num_grains;
static int num_steps;

static int* grid;
static int* next_grid;

static int final_result;

// --- Function Prototypes ---
void setup_benchmark(int argc, char *argv[]);
void run_computation();
void cleanup();

// --- Function Implementations ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s grid_size num_grains num_steps seed\n", argv[0]);
        exit(1);
    }

    grid_size  = atoi(argv[1]);
    num_grains = atoi(argv[2]);
    num_steps  = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    
    if (grid_size <= 0 || num_grains < 0 || num_steps < 0) {
        fprintf(stderr, "Invalid parameters: must be positive integers.\n");
        exit(1);
    }

    mt_seed(seed);

    size_t grid_elements = (size_t)grid_size * grid_size;
    size_t grid_bytes = grid_elements * sizeof(int);

    grid      = (int*)malloc(grid_bytes);
    next_grid = (int*)malloc(grid_bytes);

    if (!grid || !next_grid) {
        fprintf(stderr, "Memory allocation failed for grid of size %d x %d\n", grid_size, grid_size);
        exit(1);
    }
    
    memset(grid, 0, grid_bytes);

    // Randomly place grains
    for (int i = 0; i < num_grains; i++) {
        uint32_t r = mt_rand();
        int index = r % grid_elements;
        grid[index]++;
    }
}

void run_computation() {
    long long total_topples = 0;
    size_t grid_elements = (size_t)grid_size * grid_size;

    for (int step = 0; step < num_steps; step++) {
        // Clear the next_grid buffer to accumulate the new state
        memset(next_grid, 0, grid_elements * sizeof(int));
        
        for (int r = 0; r < grid_size; r++) {
            for (int c = 0; c < grid_size; c++) {
                int index = r * grid_size + c;
                int current_grains = grid[index];

                if (current_grains < 4) {
                    // Grains stay put if the pile is stable
                    next_grid[index] += current_grains;
                } else {
                    // Pile topples
                    int topples = current_grains / 4;
                    int remaining_grains = current_grains % 4;
                    total_topples += topples;
                    
                    next_grid[index] += remaining_grains;
                    
                    // Distribute toppled grains to neighbors
                    if (r > 0)               next_grid[(r - 1) * grid_size + c] += topples;
                    if (r < grid_size - 1)   next_grid[(r + 1) * grid_size + c] += topples;
                    if (c > 0)               next_grid[r * grid_size + (c - 1)] += topples;
                    if (c < grid_size - 1)   next_grid[r * grid_size + (c + 1)] += topples;
                    // Grains falling off the edge are lost
                }
            }
        }
        
        // Swap grid pointers for the next iteration
        int* temp = grid;
        grid = next_grid;
        next_grid = temp;
    }
    
    // The result is the total number of topples, a measure of the total activity.
    final_result = (int)(total_topples & 0xFFFFFFFF);
}

void cleanup() {
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

    printf("%d\n", final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
