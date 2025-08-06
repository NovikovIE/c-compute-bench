#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- Mersenne Twister (MT19937) Generator ---
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

// --- Benchmark Globals ---
int g_num_universe_elements;
int g_num_subsets;
uint64_t *g_subsets;
uint64_t g_target_universe_mask;
unsigned long long g_cover_count;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_universe_elements> <num_subsets> <seed>\n", argv[0]);
        exit(1);
    }

    g_num_universe_elements = atoi(argv[1]);
    g_num_subsets = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    if (g_num_universe_elements <= 0 || g_num_universe_elements > 64) {
        fprintf(stderr, "Error: num_universe_elements must be between 1 and 64.\n");
        exit(1);
    }
    if (g_num_subsets <= 0 || g_num_subsets > 63) { // 1ULL << 64 is UB
        fprintf(stderr, "Error: num_subsets must be between 1 and 63.\n");
        exit(1);
    }

    mt_seed(seed);

    g_subsets = (uint64_t*)malloc(g_num_subsets * sizeof(uint64_t));
    if (g_subsets == NULL) {
        fprintf(stderr, "Failed to allocate memory for subsets.\n");
        exit(1);
    }

    // Generate random subsets, represented as bitmasks
    for (int i = 0; i < g_num_subsets; i++) {
        g_subsets[i] = 0;
        // Generate a random subset with each element having a ~66.7% chance of being included
        // to make it likely that covers exist.
        for (int j = 0; j < g_num_universe_elements; j++) {
            if ((mt_rand() % 3) != 0) { 
                g_subsets[i] |= (1ULL << j);
            }
        }
    }
    
    // Create the bitmask for the complete universe set
    if (g_num_universe_elements == 64) {
        g_target_universe_mask = 0xFFFFFFFFFFFFFFFFULL;
    } else {
        g_target_universe_mask = (1ULL << g_num_universe_elements) - 1;
    }

    g_cover_count = 0;
}

void run_computation() {
    // Iterate through all 2^g_num_subsets possible sub-collections of S.
    // Each integer 'i' from 1 to 2^m-1 represents a sub-collection, where
    // the j-th bit of 'i' being set means the j-th subset is included.
    const uint64_t num_combinations = 1ULL << g_num_subsets;

    for (uint64_t i = 1; i < num_combinations; i++) {
        uint64_t current_union = 0ULL;

        // Build the union of subsets for the current combination
        for (int j = 0; j < g_num_subsets; j++) {
            if ((i >> j) & 1) {
                current_union |= g_subsets[j];
            }
        }

        // Check if the union covers the entire universe
        if (current_union == g_target_universe_mask) {
            g_cover_count++;
        }
    }
}

void cleanup() {
    free(g_subsets);
    g_subsets = NULL;
}

// --- Main Function ---
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%llu\n", g_cover_count);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
