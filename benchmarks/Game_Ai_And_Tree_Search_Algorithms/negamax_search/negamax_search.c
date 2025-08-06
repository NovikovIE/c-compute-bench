/**
 * @file negamax_search.c
 * @brief Benchmark for a Negamax tree search algorithm.
 *
 * This program simulates a core component of game AI by performing a Negamax
 * search on an implicit game tree. The tree's structure is defined by its
 * depth and branching factor. To focus on the search algorithm and memory
 * access patterns, the values of the terminal (leaf) nodes are pre-generated
 * with random scores and stored in a large array on the heap.
 *
 * The benchmark measures the time taken to perform the search itself,
 * excluding the setup time for generating the leaf node scores.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- Mersenne Twister (MT19937) PRNG --- (DO NOT MODIFY) ---
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

// --- Benchmark Globals ---
#define NEG_INFINITY -2147483647

static int SEARCH_DEPTH;
static int BRANCHING_FACTOR;

// Heap-allocated arrays for the benchmark data
static int* leaf_node_scores = NULL;
static long long* branch_multipliers = NULL;

// The final result of the computation
static int final_result;

// --- Core Algorithm ---

/**
 * @brief Recursive Negamax function.
 * @param depth The current depth in the search tree (counts down to 0).
 * @param leaf_offset The starting index in the `leaf_node_scores` array for the subtree rooted at this node.
 * @return The Negamax score for the current node.
 */
int negamax_recursive(int depth, long long leaf_offset) {
    if (depth == 0) {
        // At a leaf node, return its pre-calculated score
        return leaf_node_scores[leaf_offset];
    }

    int max_val = NEG_INFINITY;
    
    // The number of leaves covered by each child's subtree
    long long child_leaf_span = branch_multipliers[depth - 1];

    for (int i = 0; i < BRANCHING_FACTOR; ++i) {
        long long child_offset = leaf_offset + i * child_leaf_span;
        int val = -negamax_recursive(depth - 1, child_offset);
        if (val > max_val) {
            max_val = val;
        }
    }
    return max_val;
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <search_depth> <branching_factor> <seed>\n", argv[0]);
        exit(1);
    }

    SEARCH_DEPTH = atoi(argv[1]);
    BRANCHING_FACTOR = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (SEARCH_DEPTH <= 0 || BRANCHING_FACTOR <= 0) {
        fprintf(stderr, "FATAL: Search depth and branching factor must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    // Calculate total number of leaf nodes
    long long num_leaves = 1;
    for(int i = 0; i < SEARCH_DEPTH; ++i) {
        // Check for potential overflow before multiplication
        if (__builtin_mul_overflow(num_leaves, BRANCHING_FACTOR, &num_leaves)) {
            fprintf(stderr, "FATAL: Number of leaf nodes exceeds long long capacity.\n");
            exit(1);
        };
    }

    // Pre-calculate multipliers for indexing: multipliers[d] = B^d
    branch_multipliers = (long long*)malloc(SEARCH_DEPTH * sizeof(long long));
    if (branch_multipliers == NULL) {
        fprintf(stderr, "FATAL: Memory allocation for multipliers failed.\n");
        exit(1);
    }
    branch_multipliers[0] = 1;
    for (int i = 1; i < SEARCH_DEPTH; ++i) {
        branch_multipliers[i] = branch_multipliers[i - 1] * BRANCHING_FACTOR;
    }

    // Allocate and populate leaf node scores
    leaf_node_scores = (int*)malloc(num_leaves * sizeof(int));
    if (leaf_node_scores == NULL) {
        fprintf(stderr, "FATAL: Memory allocation for leaf nodes failed.\n");
        free(branch_multipliers);
        exit(1);
    }

    for (long long i = 0; i < num_leaves; ++i) {
        // Generate a score in the range [-1000, 1000]
        leaf_node_scores[i] = (int)(mt_rand() % 2001) - 1000;
    }
}

void run_computation() {
    final_result = negamax_recursive(SEARCH_DEPTH, 0);
}

void cleanup() {
    free(leaf_node_scores);
    leaf_node_scores = NULL;
    free(branch_multipliers);
    branch_multipliers = NULL;
}

// --- Main and Timing ---

int main(int argc, char* argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final result to stdout
    printf("%d\n", final_result);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
