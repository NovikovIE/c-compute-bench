#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- Global data structures ---
int program_branching_factor;
int max_path_depth;
uint32_t** constraint_table;
unsigned long long final_result;

// --- Mersenne Twister (MT19937) PRNG ---
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

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <program_branching_factor> <max_path_depth> <seed>\n", argv[0]);
        exit(1);
    }

    program_branching_factor = atoi(argv[1]);
    max_path_depth = atoi(argv[2]);
    uint32_t seed = (uint32_t)strtoul(argv[3], NULL, 10);

    if (program_branching_factor <= 0 || max_path_depth <= 0) {
        fprintf(stderr, "FATAL: Branching factor and path depth must be positive integers.\n");
        exit(1);
    }

    mt_seed(seed);
    
    // Allocate a 2D array to store pre-computed constraints for each path decision.
    constraint_table = (uint32_t**)malloc(max_path_depth * sizeof(uint32_t*));
    if (constraint_table == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed for constraint_table rows.\n");
        exit(1);
    }

    for (int i = 0; i < max_path_depth; ++i) {
        constraint_table[i] = (uint32_t*)malloc(program_branching_factor * sizeof(uint32_t));
        if (constraint_table[i] == NULL) {
            fprintf(stderr, "FATAL: Memory allocation failed for constraint_table columns.\n");
            exit(1);
        }
        for (int j = 0; j < program_branching_factor; ++j) {
            constraint_table[i][j] = mt_rand();
        }
    }
    final_result = 0;
}

// Helper function to recursively explore the state space
static unsigned long long explore_paths(int depth, unsigned long long path_value) {
    if (depth >= max_path_depth) {
        return path_value;
    }

    unsigned long long sub_tree_sum = 0;
    for (int i = 0; i < program_branching_factor; ++i) {
        // Combine current path value with a new constraint from the pre-generated table.
        // This simulates adding a new logical conjunct to the path condition.
        unsigned long long new_path_value = path_value ^ constraint_table[depth][i];
        
        // Add a simple, pseudo-random transformation to simulate constraint solver work.
        // This is a 64-bit mixing function to make the work non-trivial.
        new_path_value *= 0x9E3779B97F4A7C15ULL;
        new_path_value ^= new_path_value >> 32;

        sub_tree_sum += explore_paths(depth + 1, new_path_value);
    }
    return sub_tree_sum;
}

void run_computation() {
    // Start symbolic execution from the root (depth 0) with an initial path value.
    // The total number of paths explored is program_branching_factor ^ max_path_depth.
    final_result = explore_paths(0, 1ULL);
}

void cleanup() {
    if (constraint_table != NULL) {
        for (int i = 0; i < max_path_depth; ++i) {
            free(constraint_table[i]);
        }
        free(constraint_table);
    }
}


int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final accumulated result to stdout to prevent dead code elimination
    printf("%llu\n", final_result);

    // Print the execution time to stderr
    fprintf(stderr, "%.6f", time_taken);
    
    return 0;
}
