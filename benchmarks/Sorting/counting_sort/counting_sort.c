#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

// --- Mersenne Twister (MT19937) start ---
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
// --- Mersenne Twister (MT19937) end ---

// --- Benchmark Globals ---
int* data_to_sort;
int* sorted_data;
int* count_array;
int num_elements;
int max_value_in_data;
int final_result; // For preventing dead code elimination

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s num_elements max_value_in_data seed\n", argv[0]);
        exit(1);
    }

    num_elements = atoi(argv[1]);
    max_value_in_data = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (num_elements <= 0 || max_value_in_data <= 0) {
        fprintf(stderr, "FATAL: num_elements and max_value_in_data must be positive integers.\n");
        exit(1);
    }

    mt_seed(seed);

    data_to_sort = (int*)malloc(num_elements * sizeof(int));
    sorted_data = (int*)malloc(num_elements * sizeof(int));
    // count array needs to hold indices from 0 to max_value_in_data
    count_array = (int*)malloc((max_value_in_data + 1) * sizeof(int));

    if (!data_to_sort || !sorted_data || !count_array) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    for (int i = 0; i < num_elements; ++i) {
        data_to_sort[i] = mt_rand() % (max_value_in_data + 1);
    }
}

void run_computation() {
    // Step 1: Initialize count array to all zeros.
    memset(count_array, 0, (max_value_in_data + 1) * sizeof(int));

    // Step 2: Store the count of each unique element.
    for (int i = 0; i < num_elements; ++i) {
        count_array[data_to_sort[i]]++;
    }

    // Step 3: Store the cumulative count of each element.
    for (int i = 1; i <= max_value_in_data; ++i) {
        count_array[i] += count_array[i - 1];
    }

    // Step 4: Build the output array from the count array.
    // Iterate from the end to maintain stability (for sorting objects).
    for (int i = num_elements - 1; i >= 0; --i) {
        sorted_data[count_array[data_to_sort[i]] - 1] = data_to_sort[i];
        count_array[data_to_sort[i]]--;
    }

    // Calculate a final result to prevent dead code elimination.
    long long accumulator = 0;
    int checksum_count = num_elements > 1000 ? 1000 : num_elements;
    for (int i = 0; i < checksum_count; ++i) {
        accumulator += sorted_data[i];
    }
    final_result = (int)(accumulator & 0xFFFFFFFF);
}

void cleanup() {
    free(data_to_sort);
    free(sorted_data);
    free(count_array);
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
    printf("%d\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
