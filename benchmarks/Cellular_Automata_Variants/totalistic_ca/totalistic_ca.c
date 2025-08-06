/*
 * Benchmark: totalistic_ca
 * Description: Simulates a totalistic cellular automaton on a 2D grid.
 *              The next state of a cell is determined solely by the sum of its
 *              own state and its 8 neighbors' states (Moore neighborhood).
 *              The rules are generated randomly at setup.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- Mersenne Twister (MT19937) Generator (DO NOT MODIFY) ---
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

// Benchmark parameters and data
int grid_size;
int num_steps;
int num_states;

int *grid;       // Current state grid
int *next_grid;  // Next state grid
int *rule_table; // Maps neighbor sum to next state

unsigned long long final_result; // To prevent dead-code elimination

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <grid_size> <num_steps> <num_states> <seed>\n", argv[0]);
        exit(1);
    }

    grid_size = atoi(argv[1]);
    num_steps = atoi(argv[2]);
    num_states = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    
    mt_seed(seed);

    if (grid_size <= 0 || num_steps <= 0 || num_states <= 1) {
        fprintf(stderr, "FATAL: Invalid parameters. grid_size > 0, num_steps > 0, num_states > 1\n");
        exit(1);
    }

    long long grid_mem_size = (long long)grid_size * grid_size * sizeof(int);
    grid = (int*)malloc(grid_mem_size);
    next_grid = (int*)malloc(grid_mem_size);

    if (!grid || !next_grid) {
        fprintf(stderr, "FATAL: Memory allocation failed for grids.\n");
        exit(1);
    }

    // Initialize grid with random states
    for (int i = 0; i < grid_size * grid_size; ++i) {
        grid[i] = mt_rand() % num_states;
    }

    // For a Moore neighborhood (8 neighbors + self), there are 9 cells.
    // The maximum possible sum of states is 9 * (num_states - 1).
    int max_sum = 9 * (num_states - 1);
    rule_table = (int*)malloc((max_sum + 1) * sizeof(int));
    if (!rule_table) {
        fprintf(stderr, "FATAL: Memory allocation failed for rule table.\n");
        free(grid);
        free(next_grid);
        exit(1);
    }

    // Generate a random rule table
    for (int i = 0; i <= max_sum; ++i) {
        rule_table[i] = mt_rand() % num_states;
    }
}

void run_computation() {
    for (int step = 0; step < num_steps; ++step) {
        for (int i = 0; i < grid_size; ++i) {
            for (int j = 0; j < grid_size; ++j) {
                int sum = 0;
                // Sum the Moore neighborhood (self + 8 neighbors) with wrap-around
                for (int di = -1; di <= 1; ++di) {
                    for (int dj = -1; dj <= 1; ++dj) {
                        int ni = (i + di + grid_size) % grid_size;
                        int nj = (j + dj + grid_size) % grid_size;
                        sum += grid[ni * grid_size + nj];
                    }
                }
                next_grid[i * grid_size + j] = rule_table[sum];
            }
        }

        // Swap grids for the next iteration
        int* temp = grid;
        grid = next_grid;
        next_grid = temp;
    }

    // Calculate a final checksum to ensure computation is not optimized away
    final_result = 0;
    for (int i = 0; i < grid_size * grid_size; ++i) {
        final_result += grid[i];
    }
}

void cleanup() {
    if (grid) free(grid);
    if (next_grid) free(next_grid);
    if (rule_table) free(rule_table);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
    cleanup();

    // Print result to stdout
    printf("%llu\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
