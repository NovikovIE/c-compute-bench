#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

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

// --- Benchmark Globals ---
int num_items;
int bin_capacity;
int *items;         // Array of item sizes
int *bins;          // Array of remaining space in each bin
int final_bins_used; // The result

// Comparison function for qsort to sort in descending order
int compare_desc(const void *a, const void *b) {
    int val_a = *(const int *)a;
    int val_b = *(const int *)b;
    if (val_a < val_b) return 1;
    if (val_a > val_b) return -1;
    return 0;
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_items> <bin_capacity> <seed>\n", argv[0]);
        exit(1);
    }

    num_items = atoi(argv[1]);
    bin_capacity = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (num_items <= 0 || bin_capacity <= 0) {
        fprintf(stderr, "Error: num_items and bin_capacity must be positive integers.\n");
        exit(1);
    }

    mt_seed(seed);

    items = (int *)malloc(num_items * sizeof(int));
    if (!items) {
        fprintf(stderr, "Failed to allocate memory for items.\n");
        exit(1);
    }

    // The maximum number of bins needed is num_items
    bins = (int *)malloc(num_items * sizeof(int));
    if (!bins) {
        fprintf(stderr, "Failed to allocate memory for bins.\n");
        free(items);
        exit(1);
    }

    // Generate random item sizes between 1 and bin_capacity
    for (int i = 0; i < num_items; ++i) {
        items[i] = (mt_rand() % bin_capacity) + 1;
    }
}

void run_computation() {
    // First-Fit Decreasing algorithm
    
    // 1. Sort items in decreasing order of size
    qsort(items, num_items, sizeof(int), compare_desc);

    // 2. Place items into bins
    int current_bins_count = 0;
    for (int i = 0; i < num_items; ++i) {
        int item_size = items[i];
        int j;
        // Find the first bin that can fit the item
        for (j = 0; j < current_bins_count; ++j) {
            if (bins[j] >= item_size) {
                bins[j] -= item_size;
                break;
            }
        }

        // If no bin was found, open a new one
        if (j == current_bins_count) {
            bins[current_bins_count] = bin_capacity - item_size;
            current_bins_count++;
        }
    }

    final_bins_used = current_bins_count;
}

void cleanup() {
    free(items);
    free(bins);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout
    printf("%d\n", final_bins_used);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
