#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- BEGIN: Mersenne Twister (DO NOT MODIFY) ---
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
// --- END: Mersenne Twister ---


// --- Benchmark Globals ---
// The number of pairs of up/down steps. The path length is 2*N.
static int N;
// The final result: total number of valid Dyck paths.
static unsigned long long path_count;


// Forward declaration for the recursive helper
unsigned long long count_paths_recursive(int x, int y);


/**
 * @brief Parses command line arguments, and initializes benchmark parameters.
 * The core logic does not require any data structures, so this function
 * primarily sets the N parameter for the combinatorial problem.
 */
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_pairs> <seed>\n", argv[0]);
        exit(1);
    }

    N = atoi(argv[1]);
    uint32_t seed = (uint32_t)strtoul(argv[2], NULL, 10);
    
    mt_seed(seed); // Seed the generator, even if not directly used in computation

    if (N < 0) {
        fprintf(stderr, "Error: num_pairs must be a non-negative integer.\n");
        exit(1);
    }
    if (N > 30) {
        // Unsigned long long overflows around N=34
        fprintf(stderr, "Warning: num_pairs > 30 may lead to overflow or excessive runtime.\n");
    }
}

/**
 * @brief Recursively counts the number of Dyck paths of length 2*N.
 * A Dyck path is a sequence of +1 and -1 steps that starts and ends at 0,
 * and never goes below 0.
 * @param x The number of steps taken so far (current horizontal position).
 * @param y The current altitude (current vertical position).
 * @return The number of valid paths from the current state (x, y).
 */
unsigned long long count_paths_recursive(int x, int y) {
    // Pruning rule 1: If we are below the x-axis, the path is invalid.
    if (y < 0) {
        return 0;
    }
    
    // Pruning rule 2: If we are too high to return to y=0 in the remaining steps,
    // the path is invalid. Remaining steps = 2*N - x.
    if (y > (2 * N - x)) {
        return 0;
    }

    // Base case: If we have taken 2*N steps and ended at y=0, we found a valid path.
    if (x == 2 * N) {
        return (y == 0) ? 1 : 0;
    }

    // Recursive step: sum the paths from taking an 'up' step and a 'down' step.
    return count_paths_recursive(x + 1, y + 1) + count_paths_recursive(x + 1, y - 1);
}

/**
 * @brief Executes the core computation: counting Dyck paths.
 * This function initiates the recursive counting process.
 */
void run_computation() {
    path_count = count_paths_recursive(0, 0);
}

/**
 * @brief Frees any resources allocated in setup_benchmark.
 * In this specific benchmark, no dynamic memory is used, so this is a no-op.
 */
void cleanup() {
    // No heap memory was allocated.
}


int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final result (total number of paths) to stdout.
    printf("%llu\n", path_count);

    // Print the execution time to stderr.
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
