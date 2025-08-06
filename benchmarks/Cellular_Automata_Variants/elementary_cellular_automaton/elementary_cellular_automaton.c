#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

// --- Mersenne Twister (MT19937) START ---
// Do Not Modify
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
// --- Mersenne Twister (MT19937) END ---

// Global structure for benchmark data
typedef struct {
    int rule_number;
    int grid_width;
    int num_steps;
    unsigned char *grid;
    unsigned char *next_grid;
    unsigned char rule_map[8];
    int final_result;
} BenchmarkData;

static BenchmarkData B;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <rule_number> <grid_width> <num_steps> <seed>\n", argv[0]);
        exit(1);
    }

    B.rule_number = atoi(argv[1]);
    B.grid_width = atoi(argv[2]);
    B.num_steps = atoi(argv[3]);
    uint32_t seed = (uint32_t)strtoul(argv[4], NULL, 10);

    mt_seed(seed);

    if (B.rule_number < 0 || B.rule_number > 255) {
        fprintf(stderr, "Error: rule_number must be between 0 and 255.\n");
        exit(1);
    }

    // Pre-calculate the rule's output for each of the 8 possible neighborhoods
    for (int i = 0; i < 8; ++i) {
        B.rule_map[i] = (B.rule_number >> i) & 1;
    }

    B.grid = (unsigned char *)malloc(B.grid_width * sizeof(unsigned char));
    B.next_grid = (unsigned char *)malloc(B.grid_width * sizeof(unsigned char));
    if (B.grid == NULL || B.next_grid == NULL) {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize the grid. A common setup is a single '1' in the middle.
    memset(B.grid, 0, B.grid_width * sizeof(unsigned char));
    B.grid[B.grid_width / 2] = 1;
}

void run_computation() {
    for (int step = 0; step < B.num_steps; ++step) {
        for (int i = 0; i < B.grid_width; ++i) {
            // Get neighborhood states with periodic boundary conditions
            int left_idx = (i == 0) ? B.grid_width - 1 : i - 1;
            int right_idx = (i == B.grid_width - 1) ? 0 : i + 1;

            unsigned char left = B.grid[left_idx];
            unsigned char center = B.grid[i];
            unsigned char right = B.grid[right_idx];

            // The neighborhood forms a 3-bit index
            int rule_index = (left << 2) | (center << 1) | right;
            
            // Apply the rule
            B.next_grid[i] = B.rule_map[rule_index];
        }

        // Swap grids for the next iteration
        unsigned char *temp = B.grid;
        B.grid = B.next_grid;
        B.next_grid = temp;
    }

    // Calculate a final result to prevent dead code elimination
    int sum = 0;
    for (int i = 0; i < B.grid_width; ++i) {
        sum += B.grid[i];
    }
    B.final_result = sum;
}

void cleanup() {
    free(B.grid);
    free(B.next_grid);
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
    printf("%d\n", B.final_result);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
