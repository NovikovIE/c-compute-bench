#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <string.h>

// --- Mersenne Twister (DO NOT MODIFY) ---
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
typedef struct {
    unsigned char *input_data;
    int *output_data;
    size_t data_size;
    int alphabet_size;
} BenchmarkData;

BenchmarkData g_data;
long long g_result_accumulator = 0;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <input_size_mb> <alphabet_size> <seed>\n", argv[0]);
        exit(1);
    }

    long input_size_mb = atol(argv[1]);
    g_data.alphabet_size = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if(input_size_mb <= 0 || g_data.alphabet_size <= 1 || g_data.alphabet_size > 256) {
        fprintf(stderr, "Invalid arguments: size_mb must be > 0, alphabet_size must be in [2, 256]\n");
        exit(1);
    }

    g_data.data_size = input_size_mb * 1024 * 1024;
    mt_seed(seed);

    g_data.input_data = (unsigned char*)malloc(g_data.data_size * sizeof(unsigned char));
    if (!g_data.input_data) {
        fprintf(stderr, "Failed to allocate memory for input data.\n");
        exit(1);
    }

    g_data.output_data = (int*)malloc(g_data.data_size * sizeof(int));
    if (!g_data.output_data) {
        fprintf(stderr, "Failed to allocate memory for output data.\n");
        free(g_data.input_data);
        exit(1);
    }

    // Generate random input data based on the alphabet
    for (size_t i = 0; i < g_data.data_size; ++i) {
        g_data.input_data[i] = (unsigned char)(mt_rand() % g_data.alphabet_size);
    }
}

void run_computation() {
    // Use a stack-allocated array for the alphabet; it's small (max 256).
    unsigned char alphabet[256]; 
    long long local_accumulator = 0;

    // Initialize the alphabet list
    for (int i = 0; i < g_data.alphabet_size; ++i) {
        alphabet[i] = (unsigned char)i;
    }

    for (size_t i = 0; i < g_data.data_size; ++i) {
        unsigned char current_symbol = g_data.input_data[i];
        int found_index = 0; // Default to 0, if not found loop will correct

        // 1. Find the symbol's index in the current alphabet list
        // This linear scan is the performance bottleneck, as intended.
        for (int j = 0; j < g_data.alphabet_size; ++j) {
            if (alphabet[j] == current_symbol) {
                found_index = j;
                break;
            }
        }

        // 2. Output the index and add to accumulator
        g_data.output_data[i] = found_index;
        local_accumulator += found_index;

        // 3. Move the symbol to the front of the list
        if (found_index > 0) {
            // Use memmove for safe overlapping memory block moving
            memmove(&alphabet[1], &alphabet[0], found_index * sizeof(unsigned char));
            alphabet[0] = current_symbol;
        }
    }
    
    g_result_accumulator = local_accumulator;
}

void cleanup() {
    free(g_data.input_data);
    free(g_data.output_data);
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

    printf("%lld\n", g_result_accumulator);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
