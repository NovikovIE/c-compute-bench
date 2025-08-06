#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- START: Mersenne Twister (Do Not Modify) ---
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
// --- END: Mersenne Twister ---

// Global data structure for benchmark parameters and data
typedef struct {
    int num_elements;
    int *array;
    int final_result;
} BenchmarkData;

BenchmarkData g_data;

// Finds the index of the maximum element in arr[0..n-1]
int find_max_index(int arr[], int n) {
    int max_i = 0;
    for (int i = 1; i < n; i++) {
        if (arr[i] > arr[max_i]) {
            max_i = i;
        }
    }
    return max_i;
}

// Reverses arr[0..i]
void flip(int arr[], int i) {
    int temp, start = 0;
    while (start < i) {
        temp = arr[start];
        arr[start] = arr[i];
        arr[i] = temp;
        start++;
        i--;
    }
}

// Parses arguments, allocates memory, and generates input data.
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_elements> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_elements = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if (g_data.num_elements <= 0) {
        fprintf(stderr, "FATAL: num_elements must be a positive integer.\n");
        exit(1);
    }

    mt_seed(seed);

    g_data.array = (int *)malloc(g_data.num_elements * sizeof(int));
    if (g_data.array == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    for (int i = 0; i < g_data.num_elements; i++) {
        g_data.array[i] = (int)mt_rand();
    }
}

// Executes the core pancake sort algorithm.
void run_computation() {
    for (int curr_size = g_data.num_elements; curr_size > 1; --curr_size) {
        int max_idx = find_max_index(g_data.array, curr_size);

        if (max_idx != curr_size - 1) {
            // Move the maximum element to the front
            if (max_idx != 0) {
                flip(g_data.array, max_idx);
            }
            // Flip the entire current subarray to move the maximum element to its correct final position
            flip(g_data.array, curr_size - 1);
        }
    }

    // Calculate a checksum of the sorted array to prevent dead code elimination.
    g_data.final_result = 0;
    for (int i = 0; i < g_data.num_elements; i++) {
        g_data.final_result ^= g_data.array[i];
    }
}

// Frees all memory allocated in setup_benchmark.
void cleanup() {
    free(g_data.array);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final checksum result to stdout
    printf("%d\n", g_data.final_result);

    // Print the execution time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
