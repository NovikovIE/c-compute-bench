/**
 * @file generate_combinations.c
 * @brief Benchmark for counting combinations (n-choose-k) using a recursive approach.
 *
 * This program calculates the number of ways to choose k items from a set of n
 * distinct items, denoted as C(n, k). The computation is performed using the
 * recursive Pascal's identity: C(n, k) = C(n-1, k-1) + C(n-1, k).
 * This method is computationally intensive due to its exponential number of
 * recursive calls, making it a good benchmark for CPU-bound tasks involving
 * integer arithmetic and function call overhead.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- Mersenne Twister (MT19937) PRNG --- VERBATIM AS REQUESTED ---
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

// --- Benchmark State ---

// Global struct to hold benchmark data and parameters.
// This is the primary way data is shared between setup, computation, and cleanup.
typedef struct {
    int n_total_elements;
    int k_selection_size;
    unsigned long long result;
} BenchmarkData;

static BenchmarkData g_data;

// --- Core Computation ---

/**
 * @brief Recursively calculates the number of combinations C(n, k).
 * 
 * @param n The total number of elements.
 * @param k The number of elements to choose.
 * @return The number of combinations. Result fits in unsigned long long for the given constraints.
 */
unsigned long long count_combinations_recursive(int n, int k) {
    // Base case: invalid selection
    if (k < 0 || k > n) {
        return 0;
    }
    // Base case: choosing 0 or all elements is 1 way
    if (k == 0 || k == n) {
        return 1;
    }
    // Optimization using identity C(n, k) = C(n, n-k)
    if (k > n / 2) {
        k = n - k;
    }
    // Recursive step using Pascal's identity
    return count_combinations_recursive(n - 1, k - 1) + count_combinations_recursive(n - 1, k);
}

// --- Benchmark Functions ---

/**
 * @brief Parses arguments, initializes data.
 * This benchmark is compute-bound and does not require large memory allocations.
 * The 'data' are the parameters n and k themselves.
 */
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s n_total_elements k_selection_size seed\n", argv[0]);
        exit(1);
    }

    g_data.n_total_elements = atoi(argv[1]);
    g_data.k_selection_size = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    // Seed the PRNG (not used in this specific computation, but required by spec)
    mt_seed(seed);

    // Initialize result holder
    g_data.result = 0;

    if (g_data.n_total_elements < 0 || g_data.k_selection_size < 0 || g_data.k_selection_size > g_data.n_total_elements) {
        fprintf(stderr, "Invalid input: n and k must be non-negative, and k <= n.\n");
        exit(1);
    }
}

/**
 * @brief Runs the core computation of the benchmark.
 */
void run_computation() {
    g_data.result = count_combinations_recursive(g_data.n_total_elements, g_data.k_selection_size);
}

/**
 * @brief Frees any allocated resources.
 * For this problem, no heap memory is used, so this function is empty,
 * but it is kept for structural consistency with the benchmark template.
 */
void cleanup() {
    // No memory to free.
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
    printf("%llu\n", g_data.result);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}