#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <limits.h>

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
static int SEARCH_DEPTH;
static int BRANCHING_FACTOR;
static int* value_table;
static size_t VALUE_TABLE_SIZE;
static int final_result;

// --- Helper function for run_computation ---
// This function implements the alpha-beta algorithm on a "virtual" tree.
// Leaf node values are found by hashing the path to the leaf and looking up the
// result in the pre-computed value_table.
int alphabeta(int depth, int alpha, int beta, int is_maximizing_player, uint64_t path_hash) {
    if (depth == 0) {
        // At a leaf node, find its value in the pre-generated table
        return value_table[path_hash % VALUE_TABLE_SIZE];
    }

    if (is_maximizing_player) {
        int max_eval = INT_MIN;
        for (int i = 0; i < BRANCHING_FACTOR; ++i) {
            uint64_t new_path_hash = path_hash * 31 + i; // Simple rolling hash
            int eval = alphabeta(depth - 1, alpha, beta, 0, new_path_hash);
            if(eval > max_eval) {
                max_eval = eval;
            }
            if(eval > alpha) {
                alpha = eval;
            }
            if (beta <= alpha) {
                break; // Beta cut-off
            }
        }
        return max_eval;
    } else {
        int min_eval = INT_MAX;
        for (int i = 0; i < BRANCHING_FACTOR; ++i) {
            uint64_t new_path_hash = path_hash * 31 + i; // Simple rolling hash
            int eval = alphabeta(depth - 1, alpha, beta, 1, new_path_hash);
             if(eval < min_eval) {
                min_eval = eval;
            }
            if(eval < beta) {
                beta = eval;
            }
            if (beta <= alpha) {
                break; // Alpha cut-off
            }
        }
        return min_eval;
    }
}

// --- Benchmark Setup, Computation, and Cleanup ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <search_depth> <branching_factor> <seed>\n", argv[0]);
        exit(1);
    }

    SEARCH_DEPTH = atoi(argv[1]);
    BRANCHING_FACTOR = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    mt_seed(seed);

    // Allocate a large table of random values. The tree search will use these
    // as leaf node evaluations. This separates data generation from computation.
    VALUE_TABLE_SIZE = 1 << 24; // 16,777,216 entries, ~64MB
    value_table = (int*)malloc(VALUE_TABLE_SIZE * sizeof(int));
    if (!value_table) {
        fprintf(stderr, "Failed to allocate memory for value_table.\n");
        exit(1);
    }

    for (size_t i = 0; i < VALUE_TABLE_SIZE; ++i) {
        // Generate scores between -100 and 100
        value_table[i] = (mt_rand() % 201) - 100;
    }
}

void run_computation() {
    // Start the search from the root of the virtual tree.
    // The root is a maximizing player, at the specified search depth.
    // Path hash starts at 0.
    final_result = alphabeta(SEARCH_DEPTH, INT_MIN, INT_MAX, 1, 0UL);
}

void cleanup() {
    free(value_table);
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

    // Print the final evaluation of the root node to stdout
    printf("%d\n", final_result);

    // Print the execution time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
