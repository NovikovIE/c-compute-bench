#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <string.h>

// --- START MERSENNE TWISTER (DO NOT MODIFY) ---
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

// Global variables for benchmark data
static unsigned char *g_mtf_data;       // Input: MTF-encoded data (indices)
static unsigned char *g_output_data;    // Output: Decoded data
static unsigned char *g_alphabet;       // The dynamic alphabet list for the transform
static size_t g_input_size;             // Number of elements in g_mtf_data
static int g_alphabet_size;             // Size of the alphabet
static unsigned long long g_checksum;   // Result to prevent dead-code elimination

/**
 * @brief Parses arguments, allocates memory, and generates input data.
 * 
 * `argv[1]` (input_size_mb): The size of the input data in megabytes.
 * `argv[2]` (alphabet_size): The number of unique symbols in the data.
 * `argv[3]` (seed): The seed for the random number generator.
 */
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <input_size_mb> <alphabet_size> <seed>\n", argv[0]);
        exit(1);
    }

    size_t input_size_mb = atol(argv[1]);
    g_alphabet_size = atoi(argv[2]);
    uint32_t seed = (uint32_t)strtoul(argv[3], NULL, 10);

    if (input_size_mb <= 0) {
        fprintf(stderr, "ERROR: input_size_mb must be positive.\n");
        exit(1);
    }
    if (g_alphabet_size <= 0 || g_alphabet_size > 256) {
        fprintf(stderr, "ERROR: alphabet_size must be between 1 and 256.\n");
        exit(1);
    }

    g_input_size = input_size_mb * 1024 * 1024;
    mt_seed(seed);

    // Allocate memory
    g_mtf_data = (unsigned char *)malloc(g_input_size * sizeof(unsigned char));
    g_output_data = (unsigned char *)malloc(g_input_size * sizeof(unsigned char));
    g_alphabet = (unsigned char *)malloc(g_alphabet_size * sizeof(unsigned char));

    if (!g_mtf_data || !g_output_data || !g_alphabet) {
        fprintf(stderr, "ERROR: Memory allocation failed.\n");
        exit(1);
    }

    // Generate random MTF-encoded data (indices)
    for (size_t i = 0; i < g_input_size; ++i) {
        g_mtf_data[i] = (unsigned char)(mt_rand() % g_alphabet_size);
    }

    // Initialize the alphabet list to a known state (0, 1, 2, ...)
    for (int i = 0; i < g_alphabet_size; ++i) {
        g_alphabet[i] = (unsigned char)i;
    }
}

/**
 * @brief Executes the inverse Move-to-Front transform.
 *
 * This function reconstructs the original data stream from a sequence of MTF indices.
 * For each index, it retrieves the symbol at that position in a dynamically maintained
 * alphabet list and then moves that symbol to the front of the list.
 */
void run_computation() {
    for (size_t i = 0; i < g_input_size; ++i) {
        unsigned char mtf_index = g_mtf_data[i];

        // Get the symbol at the given index
        unsigned char symbol = g_alphabet[mtf_index];
        g_output_data[i] = symbol;

        // Move the symbol to the front of the alphabet list
        // This is the performance-critical part of the algorithm.
        // memmove is used because the source and destination regions overlap.
        memmove(&g_alphabet[1], &g_alphabet[0], mtf_index * sizeof(unsigned char));
        g_alphabet[0] = symbol;
    }

    // Calculate a checksum to prevent the compiler from optimizing away the work
    g_checksum = 0;
    for (size_t i = 0; i < g_input_size; ++i) {
        g_checksum += g_output_data[i];
    }
}

/**
 * @brief Frees all memory allocated in setup_benchmark.
 */
void cleanup() {
    free(g_mtf_data);
    free(g_output_data);
    free(g_alphabet);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final checksum to stdout
    printf("%llu\n", g_checksum);

    // Print the execution time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
