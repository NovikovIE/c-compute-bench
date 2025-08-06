#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <math.h>

// --- Start of Mersenne Twister (verbatim) ---
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
// --- End of Mersenne Twister ---

// --- Timsort Implementation ---

// We use a fixed run size for simplicity in this benchmark implementation.
// A more complex Timsort calculates this dynamically.
const int RUN = 32;

// Insertion sort is used to sort the small chunks (runs).
void insertion_sort(int arr[], int left, int right) {
    for (int i = left + 1; i <= right; i++) {
        int temp = arr[i];
        int j = i - 1;
        while (j >= left && arr[j] > temp) {
            arr[j + 1] = arr[j];
            j--;
        }
        arr[j + 1] = temp;
    }
}

// The merge function to combine two sorted runs.
void merge(int arr[], int l, int m, int r) {
    int len1 = m - l + 1, len2 = r - m;
    int* left = (int*)malloc(len1 * sizeof(int));
    int* right = (int*)malloc(len2 * sizeof(int));

    if (!left || !right) {
        fprintf(stderr, "FATAL: Memory allocation failed in merge.\n");
        exit(1);
    }

    memcpy(left, &arr[l], len1 * sizeof(int));
    memcpy(right, &arr[m + 1], len2 * sizeof(int));

    int i = 0, j = 0, k = l;

    while (i < len1 && j < len2) {
        if (left[i] <= right[j]) {
            arr[k++] = left[i++];
        } else {
            arr[k++] = right[j++];
        }
    }

    while (i < len1) {
        arr[k++] = left[i++];
    }

    while (j < len2) {
        arr[k++] = right[j++];
    }

    free(left);
    free(right);
}

// The main Timsort algorithm. It first sorts small runs using insertion sort,
// then merges them in a bottom-up fashion.
void tim_sort_impl(int arr[], int n) {
    for (int i = 0; i < n; i += RUN) {
        insertion_sort(arr, i, (i + RUN - 1 < n - 1) ? (i + RUN - 1) : (n - 1));
    }

    for (int size = RUN; size < n; size = 2 * size) {
        for (int left = 0; left < n; left += 2 * size) {
            int mid = left + size - 1;
            int right = (left + 2 * size - 1 < n - 1) ? (left + 2 * size - 1) : (n - 1);
            if (mid < right) {
                merge(arr, left, mid, right);
            }
        }
    }
}

// --- Benchmark Globals ---
typedef struct {
    int num_elements;
    int* data;
    long final_result;
} Benchmark_Data;

Benchmark_Data g_data;

// Helper comparison function for qsort in setup
int compare_int(const void* a, const void* b) {
    int int_a = *((int*)a);
    int int_b = *((int*)b);
    if (int_a == int_b) return 0;
    return (int_a < int_b) ? -1 : 1;
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_elements> <initial_data_distribution> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_elements = atoi(argv[1]);
    const char* distribution = argv[2];
    uint32_t seed = atoi(argv[3]);

    if (g_data.num_elements <= 0) {
        fprintf(stderr, "FATAL: num_elements must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    g_data.data = (int*)malloc(g_data.num_elements * sizeof(int));
    if (g_data.data == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed for data array.\n");
        exit(1);
    }

    if (strcmp(distribution, "random") == 0) {
        for (int i = 0; i < g_data.num_elements; i++) {
            g_data.data[i] = mt_rand();
        }
    } else if (strcmp(distribution, "sorted") == 0) {
        for (int i = 0; i < g_data.num_elements; i++) {
            g_data.data[i] = mt_rand();
        }
        qsort(g_data.data, g_data.num_elements, sizeof(int), compare_int);
    } else if (strcmp(distribution, "reverse_sorted") == 0) {
        for (int i = 0; i < g_data.num_elements; i++) {
            g_data.data[i] = mt_rand();
        }
        qsort(g_data.data, g_data.num_elements, sizeof(int), compare_int);
        for (int i = 0; i < g_data.num_elements / 2; i++) {
            int temp = g_data.data[i];
            g_data.data[i] = g_data.data[g_data.num_elements - 1 - i];
            g_data.data[g_data.num_elements - 1 - i] = temp;
        }
    } else if (strcmp(distribution, "nearly_sorted") == 0) {
        for (int i = 0; i < g_data.num_elements; i++) {
            g_data.data[i] = i;
        }
        int swaps = (int)fmax(1.0, g_data.num_elements * 0.02); // 2% of elements
        for (int i = 0; i < swaps; i++) {
            int idx1 = mt_rand() % g_data.num_elements;
            int idx2 = mt_rand() % g_data.num_elements;
            int temp = g_data.data[idx1];
            g_data.data[idx1] = g_data.data[idx2];
            g_data.data[idx2] = temp;
        }
    } else {
        fprintf(stderr, "FATAL: Unknown data distribution '%s'\n", distribution);
        exit(1);
    }
    
    g_data.final_result = 0;
}

void run_computation() {
    tim_sort_impl(g_data.data, g_data.num_elements);
    
    long checksum = 0;
    int stride = (g_data.num_elements > 1000) ? (g_data.num_elements / 1000) : 1;
    for (int i = 0; i < g_data.num_elements; i += stride) {
        checksum = (checksum + g_data.data[i]) % 2147483647;
    }
    g_data.final_result = checksum;
}

void cleanup() {
    if (g_data.data) {
        free(g_data.data);
        g_data.data = NULL;
    }
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

    printf("%ld\n", g_data.final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
