#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <string.h>

// --- Mersenne Twister (MT19937) --- (Do Not Modify)
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
// --- End Mersenne Twister ---

// Item structure for the knapsack problem
typedef struct {
    int weight;
    int value;
} Item;

// Global data structure to hold benchmark state
struct {
    int num_items;
    int max_weight;
    Item *items;
    int max_value_found;
} g_data;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_items> <max_weight> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_items = atoi(argv[1]);
    g_data.max_weight = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    if (g_data.num_items <= 0 || g_data.max_weight <= 0) {
        fprintf(stderr, "FATAL: num_items and max_weight must be positive.\n");
        exit(1);
    }
    
    // Prevent using too many items, as 2^N will overflow 64-bit integers and take eons
    if (g_data.num_items > 31) {
        fprintf(stderr, "FATAL: num_items > 31 is too large for this brute-force benchmark.\n");
        exit(1);
    }

    mt_seed(seed);

    g_data.items = (Item *)malloc(g_data.num_items * sizeof(Item));
    if (g_data.items == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Generate random items
    for (int i = 0; i < g_data.num_items; i++) {
        // Weights and values between 1 and 100
        g_data.items[i].weight = (mt_rand() % 100) + 1;
        g_data.items[i].value = (mt_rand() % 100) + 1;
    }

    g_data.max_value_found = 0;
}

void run_computation() {
    // The number of subsets is 2^num_items
    long long num_subsets = 1LL << g_data.num_items;
    int best_value = 0;

    // Iterate through all possible subsets using a bitmask approach
    for (long long i = 0; i < num_subsets; i++) {
        int current_weight = 0;
        int current_value = 0;

        // Check which items are in the current subset
        for (int j = 0; j < g_data.num_items; j++) {
            if ((i >> j) & 1) {
                current_weight += g_data.items[j].weight;
                current_value += g_data.items[j].value;
            }
        }

        // If the subset is valid (within weight capacity) and has a better value, update the maximum
        if (current_weight <= g_data.max_weight && current_value > best_value) {
            best_value = current_value;
        }
    }

    // Store result in the global struct to be printed later
    g_data.max_value_found = best_value;
}

void cleanup() {
    free(g_data.items);
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

    // Print the final result to stdout
    printf("%d\n", g_data.max_value_found);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
