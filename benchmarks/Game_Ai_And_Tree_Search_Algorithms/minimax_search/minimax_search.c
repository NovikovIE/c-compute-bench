#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <stdbool.h>
#include <limits.h>

// --- Mersenne Twister (MT19937) ---
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


// --- Benchmark Data and Globals ---
typedef struct {
    int search_depth;
    int branching_factor;
    long long num_leaf_nodes;
    int* leaf_scores;
} BenchmarkData;

static BenchmarkData g_data;
static int g_final_result;
static long long leaf_node_counter;

// --- Helper Functions ---
long long power(int base, int exp) {
    long long res = 1;
    for (int i = 0; i < exp; ++i) {
        // Use compiler builtin for overflow detection for safety
        if (__builtin_mul_overflow(res, base, &res)) {
            fprintf(stderr, "FATAL: integer overflow in power calculation.\n");
            exit(1);
        }
    }
    return res;
}

// --- Core Algorithm ---
int minimax(int depth, bool is_maximizing_player) {
    // If this is a leaf node (max depth reached), return its pre-computed score
    if (depth == 0) {
        if(leaf_node_counter >= g_data.num_leaf_nodes) {
            fprintf(stderr, "FATAL: Accessed out of bounds leaf node.\n");
            exit(1);
        }
        return g_data.leaf_scores[leaf_node_counter++];
    }

    if (is_maximizing_player) {
        int best_value = INT_MIN;
        for (int i = 0; i < g_data.branching_factor; ++i) {
            int value = minimax(depth - 1, false);
            if (value > best_value) {
                best_value = value;
            }
        }
        return best_value;
    } else { // Minimizing player
        int best_value = INT_MAX;
        for (int i = 0; i < g_data.branching_factor; ++i) {
            int value = minimax(depth - 1, true);
            if (value < best_value) {
                best_value = value;
            }
        }
        return best_value;
    }
}


// --- Benchmark Functions ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <search_depth> <branching_factor> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.search_depth = atoi(argv[1]);
    g_data.branching_factor = atoi(argv[2]);
    uint32_t seed = (uint32_t)strtoul(argv[3], NULL, 10);
    
    if (g_data.search_depth <= 0 || g_data.branching_factor <= 0) {
        fprintf(stderr, "FATAL: search_depth and branching_factor must be positive integers.\n");
        exit(1);
    }

    g_data.num_leaf_nodes = power(g_data.branching_factor, g_data.search_depth);
    
    g_data.leaf_scores = (int*)malloc(g_data.num_leaf_nodes * sizeof(int));
    if (g_data.leaf_scores == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed for leaf_scores.\n");
        exit(1);
    }

    mt_seed(seed);
    for (long long i = 0; i < g_data.num_leaf_nodes; ++i) {
        // Generate scores between -1000 and 1000
        g_data.leaf_scores[i] = (int)(mt_rand() % 2001) - 1000;
    }
}

void run_computation() {
    leaf_node_counter = 0;
    // Start the search from the root node (maximizing player)
    g_final_result = minimax(g_data.search_depth, true);
}

void cleanup() {
    free(g_data.leaf_scores);
    g_data.leaf_scores = NULL;
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

    // Print the final result to stdout
    printf("%d\n", g_final_result);

    // Print the execution time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
