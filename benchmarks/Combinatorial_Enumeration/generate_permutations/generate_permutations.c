#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// START of Mersenne Twister (DO NOT MODIFY)
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
// END of Mersenne Twister

// --- Benchmark-specific data structures and globals ---
typedef struct {
    int num_elements;
    int* data_array;
    unsigned long long checksum;
} BenchmarkData;

static BenchmarkData g_data;

// --- Forward declarations for benchmark functions ---
void generate_permutations_recursive(int k, int* arr);
void swap(int* a, int* b);

// --- Benchmark implementation ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_elements> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_elements = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if (g_data.num_elements <= 0 || g_data.num_elements > 15) { // 15! is huge
        fprintf(stderr, "Error: num_elements must be a positive integer <= 15.\n");
        exit(1);
    }

    mt_seed(seed);

    g_data.data_array = (int*)malloc(g_data.num_elements * sizeof(int));
    if (g_data.data_array == NULL) {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize array with 0, 1, 2, ..., N-1
    for (int i = 0; i < g_data.num_elements; ++i) {
        g_data.data_array[i] = i;
    }

    g_data.checksum = 0;
}

void cleanup() {
    free(g_data.data_array);
}

void swap(int* a, int* b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

// Heap's algorithm for generating permutations
void generate_permutations_recursive(int k, int* arr) {
    if (k == 1) {
        // A permutation is generated. Compute a value from it to prevent optimization.
        g_data.checksum += arr[0] ^ arr[g_data.num_elements - 1];
        return;
    }

    // Generate permutations with k-th element fixed
    generate_permutations_recursive(k - 1, arr);

    for (int i = 0; i < k - 1; i++) {
        // Swap based on whether k is odd or even
        if (k % 2 == 0) {
            swap(&arr[i], &arr[k - 1]);
        } else {
            swap(&arr[0], &arr[k - 1]);
        }
        generate_permutations_recursive(k - 1, arr);
    }
}

void run_computation() {
    if (g_data.num_elements > 0) {
        generate_permutations_recursive(g_data.num_elements, g_data.data_array);
    }
}

// --- Main function with timing ---

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%llu\n", g_data.checksum);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
