#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- Mersenne Twister (MT19937) Generator ---
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

// --- Benchmark Globals ---
int* g_data = NULL;
int* g_temp = NULL; // Auxiliary array for merging
size_t g_num_elements = 0;
int g_final_result = 0;

// --- Benchmark Functions ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_elements> <seed>\n", argv[0]);
        exit(1);
    }

    g_num_elements = atol(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if (g_num_elements == 0) {
        // Allow running with 0 elements as a valid, albeit trivial, case
        return;
    }

    // Seed the random number generator
    mt_seed(seed);

    // Allocate memory
    g_data = (int*)malloc(g_num_elements * sizeof(int));
    g_temp = (int*)malloc(g_num_elements * sizeof(int));
    if (g_data == NULL || g_temp == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Fill the array with random data
    for (size_t i = 0; i < g_num_elements; ++i) {
        g_data[i] = (int)mt_rand();
    }
}

void merge(int arr[], int temp[], size_t left, size_t mid, size_t right) {
    size_t i = left;
    size_t j = mid + 1;
    size_t k = left;

    // Copy both halves into the temporary array
    for(size_t l=left; l<=right; l++) {
        temp[l] = arr[l];
    }

    // Merge the temp array back into arr
    while (i <= mid && j <= right) {
        if (temp[i] <= temp[j]) {
            arr[k++] = temp[i++];
        } else {
            arr[k++] = temp[j++];
        }
    }

    // Copy the remaining elements of left subarray, if any
    while (i <= mid) {
        arr[k++] = temp[i++];
    }

    // Copy the remaining elements of right subarray, if any
    while (j <= right) {
        arr[k++] = temp[j++];
    }
}

void merge_sort_recursive(int arr[], int temp[], size_t left, size_t right) {
    if (left < right) {
        size_t mid = left + (right - left) / 2;
        merge_sort_recursive(arr, temp, left, mid);
        merge_sort_recursive(arr, temp, mid + 1, right);
        merge(arr, temp, left, mid, right);
    }
}

void run_computation() {
    if (g_num_elements > 0) {
        merge_sort_recursive(g_data, g_temp, 0, g_num_elements - 1);
    }
    
    // Calculate a result to prevent dead code elimination.
    // XOR-sum of a few elements.
    if (g_num_elements == 0) {
        g_final_result = 0;
    } else if (g_num_elements == 1) {
        g_final_result = g_data[0];
    } else if (g_num_elements == 2) {
        g_final_result = g_data[0] ^ g_data[g_num_elements-1];
    } else {
        g_final_result = g_data[0] ^ g_data[g_num_elements / 2] ^ g_data[g_num_elements - 1];
    }
}

void cleanup() {
    free(g_data);
    free(g_temp);
    g_data = NULL;
    g_temp = NULL;
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
    printf("%d\n", g_final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
