#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- BEGIN MERSENNE TWISTER ---
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
// --- END MERSENNE TWISTER ---

// --- BENCHMARK DATA ---
// Parameters
int MAX_DEPTH;
int BRANCHING_FACTOR;
// Result accumulator
long long final_result;

// --- ALGORITHM ---

// Recursive Depth-Limited Search (DLS)
// Simulates traversing a virtual game tree. Nodes are not stored in memory;
// their state is represented by their path, encoded in 'path_hash'.
// The value of a leaf node is derived from its path hash.
static long long dls_search(int current_depth, int depth_limit, uint64_t path_hash) {
    // Base case: If we've reached the depth limit, this is a leaf node.
    // Evaluate it and return its score. The evaluation is a simple function
    // of the path to ensure the work is not optimized away.
    if (current_depth == depth_limit) {
        return (path_hash % 997) + (path_hash % 199) + 1;
    }

    long long current_score = 0;

    // Recursive step: Visit all children.
    for (int i = 0; i < BRANCHING_FACTOR; ++i) {
        // Generate a unique hash for the child node based on the parent's hash and the child index.
        uint64_t child_path_hash = (path_hash * 313 + i);
        current_score += dls_search(current_depth + 1, depth_limit, child_path_hash);
    }

    return current_score;
}

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <max_depth> <branching_factor> <seed>\n", argv[0]);
        exit(1);
    }
    
    MAX_DEPTH = atoi(argv[1]);
    BRANCHING_FACTOR = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    if (MAX_DEPTH <= 0 || BRANCHING_FACTOR <= 0) {
        fprintf(stderr, "FATAL: max_depth and branching_factor must be positive integers.\n");
        exit(1);
    }
    
    // Seed the random number generator (for compliance, not used in core logic)
    mt_seed(seed);
    
    final_result = 0;

    // This benchmark uses a virtual tree, so no large data structures are
    // allocated on the heap. The tree is explored on-the-fly, which is
    // a common and memory-efficient strategy for vast search spaces.
}

void run_computation() {
    long long total_score = 0;

    // Iterative Deepening: Call Depth-Limited Search with increasing depth.
    for (int d = 1; d <= MAX_DEPTH; ++d) {
        // The root node starts at depth 0 with an initial path hash of 1.
        total_score += dls_search(0, d, 1);
    }
    final_result = total_score;
}

void cleanup() {
    // No memory was allocated on the heap in setup_benchmark, so nothing to free.
}

// --- MAIN ---
int main(int argc, char* argv[]) {
    struct timespec start, end;
    
    setup_benchmark(argc, argv);
    
    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);
    
    cleanup();
    
    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
    // Print result to stdout
    printf("%lld\n", final_result);
    
    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);
    
    return 0;
}