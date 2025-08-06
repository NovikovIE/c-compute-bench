#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- Mersenne Twister (MT19937) --- Do Not Modify ---
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
int* array_to_sort;
int num_elements;
int cutoff_size;
int final_result;

// --- Helper Functions ---

void swap(int* a, int* b) {
    int t = *a;
    *a = *b;
    *b = t;
}

void insertion_sort(int* arr, int low, int high) {
    for (int i = low + 1; i <= high; i++) {
        int key = arr[i];
        int j = i - 1;
        while (j >= low && arr[j] > key) {
            arr[j + 1] = arr[j];
            j = j - 1;
        }
        arr[j + 1] = key;
    }
}

// Lomuto partition scheme
int partition(int* arr, int low, int high) {
    int pivot = arr[high];
    int i = (low - 1);

    for (int j = low; j <= high - 1; j++) {
        if (arr[j] < pivot) {
            i++;
            swap(&arr[i], &arr[j]);
        }
    }
    swap(&arr[i + 1], &arr[high]);
    return (i + 1);
}

void quicksort_recursive(int* arr, int low, int high) {
    if (low < high) {
        if ((high - low + 1) < cutoff_size) {
            insertion_sort(arr, low, high);
        } else {
            int pi = partition(arr, low, high);
            quicksort_recursive(arr, low, pi - 1);
            quicksort_recursive(arr, pi + 1, high);
        }
    }
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_elements> <cutoff_size> <seed>\n", argv[0]);
        exit(1);
    }

    num_elements = atoi(argv[1]);
    cutoff_size = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (num_elements <= 0 || cutoff_size <= 0) {
        fprintf(stderr, "Invalid arguments: num_elements and cutoff_size must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    array_to_sort = (int*)malloc(num_elements * sizeof(int));
    if (array_to_sort == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(1);
    }

    for (int i = 0; i < num_elements; i++) {
        array_to_sort[i] = mt_rand();
    }
}

void run_computation() {
    quicksort_recursive(array_to_sort, 0, num_elements - 1);

    // Calculate a checksum to prevent dead code elimination.
    final_result = 0;
    // Stride to prevent this loop from dominating the benchmark time.
    for (int i = 0; i < num_elements; i += 16) {
        final_result ^= array_to_sort[i];
    }
}

void cleanup() {
    free(array_to_sort);
}

// --- Main Execution ---

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final computational result to stdout
    printf("%d\n", final_result);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
