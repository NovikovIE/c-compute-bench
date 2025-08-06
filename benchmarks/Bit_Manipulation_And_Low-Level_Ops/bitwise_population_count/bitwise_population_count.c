#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

//
// bitwise_population_count: Bit Manipulation and Low-Level Ops Benchmark
// Description: This benchmark calculates the population count (number of set bits)
// for a large array of random integers. The core computation uses the
// Brian Kernighan's algorithm, a classic bit manipulation technique.
//

// --- Mersenne Twister (MT19937) Generator (Verbatim) ---
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

// --- Global Benchmark Data ---

// Parameters
long g_num_integers;
int g_bit_width;

// Data buffer
uint32_t *g_data_array;

// Result accumulator
unsigned long long g_total_pop_count;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_integers> <bit_width> <seed>\n", argv[0]);
        exit(1);
    }

    g_num_integers = atol(argv[1]);
    g_bit_width = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (g_num_integers <= 0 || g_bit_width <= 0 || g_bit_width > 32) {
        fprintf(stderr, "Invalid arguments: num_integers must be > 0, bit_width must be between 1 and 32.\n");
        exit(1);
    }

    // Seed the random number generator
    mt_seed(seed);

    // Allocate memory on the heap
    g_data_array = (uint32_t *)malloc(g_num_integers * sizeof(uint32_t));
    if (g_data_array == NULL) {
        fprintf(stderr, "Failed to allocate memory for data array.\n");
        exit(1);
    }

    // Create a mask to constrain random numbers to the specified bit width
    uint32_t mask = (g_bit_width == 32) ? 0xFFFFFFFFUL : (1UL << g_bit_width) - 1;

    // Populate the array with random integers
    for (long i = 0; i < g_num_integers; ++i) {
        g_data_array[i] = mt_rand() & mask;
    }
}

void run_computation() {
    g_total_pop_count = 0;
    for (long i = 0; i < g_num_integers; ++i) {
        uint32_t n = g_data_array[i];
        // Brian Kernighan's algorithm to count set bits
        while (n > 0) {
            n &= (n - 1); // Clear the least significant bit set
            g_total_pop_count++;
        }
    }
}

void cleanup() {
    free(g_data_array);
    g_data_array = NULL;
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

    // Print the final accumulated result to stdout
    printf("%llu\n", g_total_pop_count);

    // Print the execution time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
