#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

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

// --- Benchmark Specific Code ---

// Struct to hold item properties
typedef struct {
    int value;
    int weight;
} Item;

// Global struct to hold all benchmark data
struct {
    int num_items;
    int capacity;
    Item *items;
    long *dp_table; // Space-optimized DP table O(capacity)
    long result;
} g_data;

// Inline helper function for max
static inline long max(long a, long b) {
    return a > b ? a : b;
}

// Setup: parse arguments, allocate memory, generate data
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_items> <capacity> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_items = atoi(argv[1]);
    g_data.capacity = atoi(argv[2]);
    uint32_t seed = (uint32_t)atol(argv[3]);

    mt_seed(seed);

    // Allocate memory for items
    g_data.items = (Item *)malloc(g_data.num_items * sizeof(Item));
    if (g_data.items == NULL) {
        fprintf(stderr, "Error: Failed to allocate memory for items.\n");
        exit(1);
    }

    // Allocate memory for the DP table (space-optimized)
    g_data.dp_table = (long *)calloc(g_data.capacity + 1, sizeof(long));
    if (g_data.dp_table == NULL) {
        fprintf(stderr, "Error: Failed to allocate memory for DP table.\n");
        free(g_data.items);
        exit(1);
    }

    // Generate random items
    for (int i = 0; i < g_data.num_items; ++i) {
        g_data.items[i].value = (mt_rand() % 1000) + 1; // value between 1 and 1000
        g_data.items[i].weight = (mt_rand() % 100) + 1; // weight between 1 and 100
    }

    g_data.result = 0;
}

// Computation: perform the 0/1 knapsack dynamic programming algorithm
void run_computation() {
    for (int i = 0; i < g_data.num_items; ++i) {
        int current_weight = g_data.items[i].weight;
        int current_value = g_data.items[i].value;
        // Iterate backwards to prevent using the same item multiple times in one pass
        for (int w = g_data.capacity; w >= current_weight; --w) {
            g_data.dp_table[w] = max(g_data.dp_table[w], current_value + g_data.dp_table[w - current_weight]);
        }
    }
    g_data.result = g_data.dp_table[g_data.capacity];
}

// Cleanup: free all allocated memory
void cleanup() {
    free(g_data.items);
    free(g_data.dp_table);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final accumulated result to stdout
    printf("%ld\n", g_data.result);

    cleanup();

    // Print timing information to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
