#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// START of Mersenne Twister 19937
// Do not modify this block
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
// END of Mersenne Twister 19937

// --- Benchmark Globals ---
typedef struct {
    int sequence_length;
    int* sequence;
    int* dp_table; 
    int final_result;
} benchmark_state_t;

static benchmark_state_t state;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <sequence_length> <seed>\n", argv[0]);
        exit(1);
    }

    state.sequence_length = atoi(argv[1]);
    uint32_t seed = atoi(argv[2]);

    if (state.sequence_length <= 0) {
        fprintf(stderr, "FATAL: sequence_length must be a positive integer.\n");
        exit(1);
    }

    mt_seed(seed);

    state.sequence = (int*)malloc(state.sequence_length * sizeof(int));
    state.dp_table = (int*)malloc(state.sequence_length * sizeof(int));

    if (!state.sequence || !state.dp_table) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    for (int i = 0; i < state.sequence_length; i++) {
        state.sequence[i] = mt_rand();
    }
}

void run_computation() {
    int n = state.sequence_length;
    if (n == 0) {
        state.final_result = 0;
        return;
    }

    // O(N^2) Dynamic Programming solution for Longest Increasing Subsequence (LIS).
    // dp_table[i] stores the length of the LIS ending at index i.
    for (int i = 0; i < n; i++) {
        state.dp_table[i] = 1; // Base case: LIS of a single element is 1.
        for (int j = 0; j < i; j++) {
            if (state.sequence[i] > state.sequence[j] && state.dp_table[i] < state.dp_table[j] + 1) {
                state.dp_table[i] = state.dp_table[j] + 1;
            }
        }
    }

    // The overall LIS length is the maximum value in the dp_table.
    int max_length = 1;
    for (int i = 0; i < n; i++) {
        if (max_length < state.dp_table[i]) {
            max_length = state.dp_table[i];
        }
    }
    state.final_result = max_length;
}

void cleanup() {
    free(state.sequence);
    free(state.dp_table);
}

// --- Main Function (with timing) ---

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final result to stdout
    printf("%d\n", state.final_result);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
