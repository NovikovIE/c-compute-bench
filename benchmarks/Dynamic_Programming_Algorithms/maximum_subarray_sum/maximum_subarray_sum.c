#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- Mersenne Twister (MT19937) Generator --- (DO NOT MODIFY)
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

// Global variables for benchmark data
static int* g_array = NULL;
static int g_array_length = 0;
static long long g_max_sum = 0;

// Function to set up the data for the benchmark
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <array_length> <seed>\n", argv[0]);
        exit(1);
    }

    g_array_length = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if (g_array_length <= 0) {
        fprintf(stderr, "FATAL: array_length must be a positive integer.\n");
        exit(1);
    }

    // Seed the random number generator
    mt_seed(seed);

    // Allocate memory for the array
    g_array = (int*)malloc(g_array_length * sizeof(int));
    if (g_array == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed for g_array.\n");
        exit(1);
    }

    // Populate the array with random integers. A mix of positive and negative
    // numbers is required for this problem to be meaningful.
    for (int i = 0; i < g_array_length; i++) {
        // Generate numbers roughly in the range [-500, 499]
        g_array[i] = (mt_rand() % 1000) - 500;
    }
}

// Function to run the core computation (Kadane's Algorithm)
void run_computation() {
    if (g_array_length == 0) {
        g_max_sum = 0;
        return;
    }

    // Kadane's algorithm is a dynamic programming approach to find the maximum subarray sum.
    // 'current_max' stores the maximum sum of a subarray ending at the current position.
    // 'global_max' stores the maximum sum found so far across the entire array.
    long long global_max = g_array[0];
    long long current_max = g_array[0];

    for (int i = 1; i < g_array_length; i++) {
        long long current_val = g_array[i];

        // The max subarray ending at `i` is either just the element at `i`, 
        // or the element at `i` combined with the max subarray ending at `i-1`.
        long long next_sum = current_max + current_val;
        current_max = (current_val > next_sum) ? current_val : next_sum;

        // Update the overall maximum sum if the new `current_max` is greater.
        if (current_max > global_max) {
            global_max = current_max;
        }
    }

    g_max_sum = global_max;
}

// Function to clean up allocated memory
void cleanup() {
    if (g_array != NULL) {
        free(g_array);
        g_array = NULL;
    }
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    cleanup();

    // Print the final result to stdout
    printf("%lld\n", g_max_sum);

    // Print the execution time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
