#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- Mersenne Twister (MT19937) Generator - DO NOT MODIFY ---
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
// --- End of Mersenne Twister ---

// Global variables for benchmark data
static int* data_array;
static size_t num_elements;
static int final_result;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_elements> <seed>\n", argv[0]);
        exit(1);
    }

    num_elements = (size_t)atol(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if (num_elements <= 0) {
        fprintf(stderr, "Error: num_elements must be a positive integer.\n");
        exit(1);
    }

    mt_seed(seed);

    data_array = (int*)malloc(num_elements * sizeof(int));
    if (data_array == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    for (size_t i = 0; i < num_elements; ++i) {
        data_array[i] = (int)mt_rand();
    }
}

void run_computation() {
    for (size_t i = 0; i < num_elements - 1; i++) {
        size_t min_idx = i;
        for (size_t j = i + 1; j < num_elements; j++) {
            if (data_array[j] < data_array[min_idx]) {
                min_idx = j;
            }
        }
        if (min_idx != i) {
            int temp = data_array[min_idx];
            data_array[min_idx] = data_array[i];
            data_array[i] = temp;
        }
    }

    // Calculate a checksum to prevent dead code elimination.
    int checksum = 0;
    for(size_t i = 0; i < num_elements; i+=16) { // Stride to reduce checksum overhead
        checksum ^= data_array[i];
    }
    final_result = checksum;
}

void cleanup() {
    free(data_array);
    data_array = NULL;
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final result to stdout to ensure computation is not optimized away
    printf("%d\n", final_result);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
