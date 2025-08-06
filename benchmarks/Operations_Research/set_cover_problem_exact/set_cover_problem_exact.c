#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- Start of Mersenne Twister (Do Not Modify) ---
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

// Global benchmark data structure
typedef struct {
    int num_universe_elements;
    int num_subsets;
    uint32_t* subsets;
    uint32_t universe_mask;
    int min_cover_size;
} BenchmarkData;

BenchmarkData g_data;

// Recursive helper function for the exact solver
void solve_recursive(int current_size, uint32_t covered_mask, int subset_idx) {
    // Pruning: if the current path is already worse than the best known solution, stop.
    if (current_size >= g_data.min_cover_size) {
        return;
    }

    // Success condition: all elements are covered. Found a new, smaller cover.
    if (covered_mask == g_data.universe_mask) {
        g_data.min_cover_size = current_size;
        return;
    }

    // Failure condition: no more subsets to try, but universe is not fully covered.
    if (subset_idx >= g_data.num_subsets) {
        return;
    }

    // --- Recursive branching ---
    
    // Branch 1: Explore solution space that INCLUDES the current subset.
    solve_recursive(current_size + 1, covered_mask | g_data.subsets[subset_idx], subset_idx + 1);

    // Branch 2: Explore solution space that EXCLUDES the current subset.
    solve_recursive(current_size, covered_mask, subset_idx + 1);
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_universe_elements> <num_subsets> <seed>\n", argv[0]);
        exit(1);
    }
    
    g_data.num_universe_elements = atoi(argv[1]);
    g_data.num_subsets = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (g_data.num_universe_elements <= 0 || g_data.num_universe_elements > 32) {
        fprintf(stderr, "FATAL: num_universe_elements must be between 1 and 32 for this bitmask implementation.\n");
        exit(1);
    }
    if (g_data.num_subsets <= 0) {
        fprintf(stderr, "FATAL: num_subsets must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    // Create the bitmask for the full universe
    g_data.universe_mask = (g_data.num_universe_elements == 32) ? 0xFFFFFFFFU : (1U << g_data.num_universe_elements) - 1;

    g_data.subsets = (uint32_t*)malloc(g_data.num_subsets * sizeof(uint32_t));
    if (!g_data.subsets) {
        fprintf(stderr, "FATAL: Memory allocation failed for subsets.\n");
        exit(1);
    }
    
    // Step 1: Generate initial random subsets to get a sparse distribution
    for (int i = 0; i < g_data.num_subsets; ++i) {
        g_data.subsets[i] = mt_rand() & mt_rand();
        g_data.subsets[i] &= g_data.universe_mask; // Keep only relevant bits
    }
    
    // Step 2: Ensure every universe element is in at least one subset to guarantee a solution exists
    uint32_t all_covered_elements = 0;
    for (int i = 0; i < g_data.num_subsets; ++i) {
        all_covered_elements |= g_data.subsets[i];
    }
    
    uint32_t uncovered_mask = g_data.universe_mask & (~all_covered_elements);
    for (int i = 0; i < g_data.num_universe_elements; ++i) {
        if ((uncovered_mask >> i) & 1) { // Check if the i-th element is uncovered
            int random_subset_idx = mt_rand() % g_data.num_subsets;
            g_data.subsets[random_subset_idx] |= (1U << i);
        }
    }
}

void run_computation() {
    // Initialize min_cover_size to a value larger than any possible solution
    g_data.min_cover_size = g_data.num_subsets + 1;
    // Start the recursive search for the exact set cover
    solve_recursive(0, 0, 0);
}

void cleanup() {
    if (g_data.subsets) {
        free(g_data.subsets);
        g_data.subsets = NULL;
    }
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // The final result is the size of the smallest cover found
    int final_result = g_data.min_cover_size;

    cleanup();

    // Print result to stdout
    printf("%d\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
