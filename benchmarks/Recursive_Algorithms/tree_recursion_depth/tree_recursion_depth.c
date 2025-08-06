#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- Mersenne Twister (verbatim) ---
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
// --- End Mersenne Twister ---


// --- Benchmark Globals ---
int PARAM_DEPTH;
long long computed_result;
// A dummy pointer to satisfy the heap allocation requirement.
char *dummy_data;
#define DUMMY_SIZE (1 * 1024 * 1024) // 1MB


// --- Core Recursive Function ---
// This function computes a sum over a ternary tree of a given depth.
// It's designed to heavily stress the function call stack.
long long tree_recursion(int current_depth, int salt) {
    if (current_depth <= 0) {
        return 1;
    }

    // Add some simple, path-dependent work at each node to prevent over-optimization.
    long long sum = (long long)((current_depth + salt) % 7);

    // Three recursive calls to create a ternary tree structure.
    sum += tree_recursion(current_depth - 1, salt + 1);
    sum += tree_recursion(current_depth - 1, salt + 2);
    sum += tree_recursion(current_depth - 1, salt + 3);

    return sum;
}


// --- Benchmark Functions ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <depth> <seed>\n", argv[0]);
        exit(1);
    }

    PARAM_DEPTH = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);
    
    mt_seed(seed);

    // Allocate and populate dummy data to ensure setup phase has work,
    // satisfying the requirement for heap allocation.
    dummy_data = (char*)malloc(DUMMY_SIZE);
    if (dummy_data == NULL) {
        fprintf(stderr, "FATAL: Failed to allocate dummy memory.\n");
        exit(1);
    }
    for (size_t i = 0; i < DUMMY_SIZE; ++i) {
        dummy_data[i] = mt_rand() & 0xFF;
    }
    
    computed_result = 0;
}

void run_computation() {
    computed_result = tree_recursion(PARAM_DEPTH, 0);
}

void cleanup() {
    free(dummy_data);
    dummy_data = NULL;
}

// --- Main --- 
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%lld\n", computed_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
