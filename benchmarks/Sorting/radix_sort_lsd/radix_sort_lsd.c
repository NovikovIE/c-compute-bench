#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <string.h>

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

// Global variables for benchmark data
static unsigned int *data_array;
static unsigned int *output_array; // Auxiliary buffer for sorting
static int num_elements;
static unsigned int max_value_in_data;
static int final_result;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_elements> <max_value_in_data> <seed>\n", argv[0]);
        exit(1);
    }

    num_elements = atoi(argv[1]);
    max_value_in_data = (unsigned int)atol(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (num_elements <= 0) {
        fprintf(stderr, "FATAL: num_elements must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    data_array = (unsigned int *)malloc(num_elements * sizeof(unsigned int));
    output_array = (unsigned int *)malloc(num_elements * sizeof(unsigned int));
    if (data_array == NULL || output_array == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    for (int i = 0; i < num_elements; ++i) {
        data_array[i] = mt_rand() % (max_value_in_data + 1);
    }
}

// Helper for radix sort: performs a stable counting sort for a given digit (byte).
static void counting_sort_pass(int shift) {
    int i;
    int count[256] = {0};

    // Store count of occurrences in count[]
    for (i = 0; i < num_elements; i++) {
        count[(data_array[i] >> shift) & 0xFF]++;
    }

    // Change count[i] so that count[i] now contains actual
    // position of this digit in output_array[]
    for (i = 1; i < 256; i++) {
        count[i] += count[i - 1];
    }

    // Build the output array
    for (i = num_elements - 1; i >= 0; i--) {
        output_array[count[(data_array[i] >> shift) & 0xFF] - 1] = data_array[i];
        count[(data_array[i] >> shift) & 0xFF]--;
    }

    // Copy the output array to data_array[], so that data_array[] now
    // contains sorted numbers according to the current digit
    memcpy(data_array, output_array, num_elements * sizeof(unsigned int));
}

void run_computation() {
    // LSD Radix Sort for 32-bit unsigned integers.
    // We process the number byte by byte, which means 4 passes for 4 bytes.
    for (int shift = 0; shift < 32; shift += 8) {
        counting_sort_pass(shift);
    }

    // To prevent dead code elimination, calculate a checksum of the sorted array.
    // A simple XOR sum of elements at a fixed stride will suffice.
    unsigned int checksum = 0;
    int stride = (num_elements > 100) ? num_elements / 100 : 1;
    for (int i = 0; i < num_elements; i += stride) {
        checksum ^= data_array[i];
    }
    final_result = (int)checksum;
}

void cleanup() {
    free(data_array);
    free(output_array);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final result to stdout to ensure computation is not optimized away.
    printf("%d\n", final_result);

    // Print the timing information to stderr.
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
