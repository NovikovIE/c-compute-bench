#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- Mersenne Twister (MT19937) Generator (Verbatim) ---
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
typedef struct {
    int grid_size;
    int block_size;
    int num_steps;
    uint8_t *grid;      // Current state of the automaton
    uint8_t *next_grid; // Buffer for the next state
    long final_result;  // Sum of all cells in the final grid
} benchmark_data_t;

static benchmark_data_t g_data;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s grid_size block_size num_steps seed\n", argv[0]);
        exit(1);
    }

    g_data.grid_size = atoi(argv[1]);
    g_data.block_size = atoi(argv[2]);
    g_data.num_steps = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    if (g_data.grid_size <= 0 || g_data.block_size <= 0 || g_data.num_steps < 0) {
        fprintf(stderr, "FATAL: grid_size, block_size, and num_steps must be positive.\n");
        exit(1);
    }
    if (g_data.block_size % 2 == 0) {
        fprintf(stderr, "FATAL: block_size must be odd for a symmetric neighborhood.\n");
        exit(1);
    }

    mt_seed(seed);

    size_t grid_bytes = (size_t)g_data.grid_size * g_data.grid_size * sizeof(uint8_t);
    g_data.grid = (uint8_t *)malloc(grid_bytes);
    g_data.next_grid = (uint8_t *)malloc(grid_bytes);

    if (!g_data.grid || !g_data.next_grid) {
        fprintf(stderr, "FATAL: Memory allocation failed for grids.\n");
        exit(1);
    }

    // Initialize the grid with random 0s and 1s
    for (int i = 0; i < g_data.grid_size * g_data.grid_size; ++i) {
        g_data.grid[i] = mt_rand() % 2;
    }

    g_data.final_result = 0;
}

void run_computation() {
    int grid_size = g_data.grid_size;
    int block_size = g_data.block_size;
    int block_offset = block_size / 2;

    // Main simulation loop
    for (int step = 0; step < g_data.num_steps; ++step) {
        // Iterate over each cell in the grid
        for (int i = 0; i < grid_size; ++i) {
            for (int j = 0; j < grid_size; ++j) {
                int sum = 0;
                // Iterate over the local block/neighborhood
                for (int bi = 0; bi < block_size; ++bi) {
                    for (int bj = 0; bj < block_size; ++bj) {
                        // Calculate neighbor coordinates with wrapping (toroidal grid)
                        int ni = (i + bi - block_offset + grid_size) % grid_size;
                        int nj = (j + bj - block_offset + grid_size) % grid_size;
                        sum += g_data.grid[ni * grid_size + nj];
                    }
                }
                // Cellular Automaton Rule: "Parity Rule"
                // The next state of a cell is 1 if the sum of states in its
                // neighborhood is odd, and 0 otherwise.
                g_data.next_grid[i * grid_size + j] = (sum % 2);
            }
        }

        // Swap grids for the next iteration
        uint8_t *temp = g_data.grid;
        g_data.grid = g_data.next_grid;
        g_data.next_grid = temp;
    }

    // Calculate a final result to prevent dead code elimination
    long total_sum = 0;
    for (int i = 0; i < grid_size * grid_size; ++i) {
        total_sum += g_data.grid[i];
    }
    g_data.final_result = total_sum;
}

void cleanup() {
    free(g_data.grid);
    free(g_data.next_grid);
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
    printf("%ld\n", g_data.final_result);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
