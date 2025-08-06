#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

// --- MERSENNE TWISTER (DO NOT MODIFY) ---
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

// Benchmark state
typedef struct {
    int block_size; // in bytes
    unsigned char *input_data;
    unsigned char *input_data_doubled;
    unsigned char *output_bwt;
    int *suffix_array;
    int final_result;
} BenchmarkData;

BenchmarkData g_data;

// Comparison function for qsort to sort cyclic shifts.
// It achieves this by comparing substrings in a doubled version of the input string,
// which elegantly handles the cyclic nature of the rotations.
int compare_rotations(const void *a, const void *b) {
    int idx1 = *(const int *)a;
    int idx2 = *(const int *)b;
    return memcmp(g_data.input_data_doubled + idx1, g_data.input_data_doubled + idx2, g_data.block_size);
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <block_size_kb> <seed>\n", argv[0]);
        exit(1);
    }

    int block_size_kb = atoi(argv[1]);
    unsigned int seed = atoi(argv[2]);

    g_data.block_size = block_size_kb * 1024;
    if (g_data.block_size <= 0) {
        fprintf(stderr, "Error: Invalid block size.\n");
        exit(1);
    }
    
    mt_seed(seed);

    // Allocate memory
    g_data.input_data = (unsigned char *)malloc(g_data.block_size * sizeof(unsigned char));
    g_data.input_data_doubled = (unsigned char *)malloc(2 * g_data.block_size * sizeof(unsigned char));
    g_data.output_bwt = (unsigned char *)malloc(g_data.block_size * sizeof(unsigned char));
    g_data.suffix_array = (int *)malloc(g_data.block_size * sizeof(int));
    g_data.final_result = 0;

    if (!g_data.input_data || !g_data.input_data_doubled || !g_data.output_bwt || !g_data.suffix_array) {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        free(g_data.input_data);
        free(g_data.input_data_doubled);
        free(g_data.output_bwt);
        free(g_data.suffix_array);
        exit(1);
    }

    // Generate random input data
    for (int i = 0; i < g_data.block_size; i++) {
        g_data.input_data[i] = mt_rand() % 256;
    }

    // Create the doubled input string to simplify cyclic comparisons
    memcpy(g_data.input_data_doubled, g_data.input_data, g_data.block_size);
    memcpy(g_data.input_data_doubled + g_data.block_size, g_data.input_data, g_data.block_size);
}

void run_computation() {
    // Initialize suffix_array with indices 0, 1, ..., N-1
    for (int i = 0; i < g_data.block_size; i++) {
        g_data.suffix_array[i] = i;
    }

    // Sort the suffixes (cyclic shifts) using qsort.
    // This is the main computational workload.
    qsort(g_data.suffix_array, g_data.block_size, sizeof(int), compare_rotations);

    // Find the original string's index in the sorted list of rotations.
    int original_index = -1;

    // Construct the BWT output string (last column of the sorted matrix)
    // and find the original index in one pass.
    for (int i = 0; i < g_data.block_size; i++) {
        int suffix_start = g_data.suffix_array[i];
        if (suffix_start == 0) {
            original_index = i;
        }
        g_data.output_bwt[i] = g_data.input_data[(suffix_start + g_data.block_size - 1) % g_data.block_size];
    }

    // Calculate a checksum of the output to prevent dead code elimination
    // and combine it with the original_index for a final result.
    long checksum = 0;
    for (int i = 0; i < g_data.block_size; i++) {
        checksum += g_data.output_bwt[i];
    }

    g_data.final_result = original_index + (int)(checksum & 0xFFFF);
}

void cleanup() {
    free(g_data.input_data);
    free(g_data.input_data_doubled);
    free(g_data.output_bwt);
    free(g_data.suffix_array);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%d\n", g_data.final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
