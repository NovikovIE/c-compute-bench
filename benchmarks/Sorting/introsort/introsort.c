#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <string.h>

// --- Mersenne Twister (MT19937) ---
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
// --- End of MT19937 ---

// --- Benchmark Globals ---
typedef struct {
    int *data;
    long num_elements;
    int max_recursion_depth;
    int final_result; // To prevent dead code elimination
} BenchmarkData;

static BenchmarkData g_data;
const int INSERTION_SORT_THRESHOLD = 16;

// --- Forward declarations for introsort helpers ---
void introsort_internal(int *begin, int *end, int depth_limit);
void insertion_sort(int *begin, int *end);
void heap_sort(int *begin, int *end);

// --- Benchmark Functions ---

void swap(int *a, int *b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

// Pivot selection for Quicksort
int* median_of_three(int *a, int *b, int *c) {
    if (*a < *b) {
        if (*b < *c) return b; // a < b < c
        if (*a < *c) return c; // a < c <= b
        return a;             // c <= a < b
    } else { // a >= b
        if (*a < *c) return a; // b <= a < c
        if (*b < *c) return c; // b < c <= a
        return b;             // c <= b <= a
    }
}

// Partition for Quicksort
int* partition_qs(int *begin, int *end) {
    int *pivot_ptr = median_of_three(begin, begin + (end - begin) / 2, end - 1);
    int pivot_value = *pivot_ptr;
    swap(pivot_ptr, end - 1); // Move pivot to end

    int *i = begin;
    for (int *j = begin; j < end - 1; j++) {
        if (*j < pivot_value) {
            swap(i, j);
            i++;
        }
    }
    swap(i, end - 1); // Move pivot to its final place
    return i;
}

// Insertion Sort for small partitions
void insertion_sort(int *begin, int *end) {
    for (int *i = begin + 1; i < end; i++) {
        int key = *i;
        int *j = i - 1;
        while (j >= begin && *j > key) {
            *(j + 1) = *j;
            j--;
        }
        *(j + 1) = key;
    }
}

// Heapify a subtree with root at index i
void heapify(int arr[], int n, int i) {
    int largest = i;
    int left = 2 * i + 1;
    int right = 2 * i + 2;

    if (left < n && arr[left] > arr[largest])
        largest = left;
    if (right < n && arr[right] > arr[largest])
        largest = right;

    if (largest != i) {
        swap(&arr[i], &arr[largest]);
        heapify(arr, n, largest);
    }
}

// Heapsort implementation
void heap_sort(int *begin, int *end) {
    int n = end - begin;
    // Build heap (rearrange array)
    for (int i = n / 2 - 1; i >= 0; i--)
        heapify(begin, n, i);

    // One by one extract an element from a heap
    for (int i = n - 1; i > 0; i--) {
        swap(&begin[0], &begin[i]);
        heapify(begin, i, 0);
    }
}

// Main Introsort recursive function
void introsort_internal(int *begin, int *end, int depth_limit) {
    while (end - begin > INSERTION_SORT_THRESHOLD) {
        if (depth_limit == 0) {
            heap_sort(begin, end);
            return;
        }
        depth_limit--;
        int *pivot = partition_qs(begin, end);
        introsort_internal(pivot + 1, end, depth_limit);
        end = pivot;
    }
    insertion_sort(begin, end);
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_elements> <max_recursion_depth> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_elements = atol(argv[1]);
    g_data.max_recursion_depth = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (g_data.num_elements <= 0) {
         fprintf(stderr, "FATAL: num_elements must be a positive integer.\n");
         exit(1);
    }
    if (g_data.max_recursion_depth <= 0) {
         fprintf(stderr, "FATAL: max_recursion_depth must be a positive integer.\n");
         exit(1);
    }

    mt_seed(seed);

    g_data.data = (int*)malloc(g_data.num_elements * sizeof(int));
    if (g_data.data == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    for (long i = 0; i < g_data.num_elements; i++) {
        g_data.data[i] = (int)mt_rand();
    }
    g_data.final_result = 0;
}

void run_computation() {
    introsort_internal(g_data.data, g_data.data + g_data.num_elements, g_data.max_recursion_depth);

    // Calculate a checksum to prevent dead code elimination
    int checksum = 0;
    long step = g_data.num_elements > 1000 ? g_data.num_elements / 1000 : 1;
    for (long i = 0; i < g_data.num_elements; i += step) {
        checksum ^= g_data.data[i];
    }
    g_data.final_result = checksum;
}

void cleanup() {
    free(g_data.data);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%d\n", g_data.final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
