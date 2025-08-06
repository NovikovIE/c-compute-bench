#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- MERSENNE TWISTER (DO NOT MODIFY) ---
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
// --- END MERSENNE TWISTER ---

// --- BENCHMARK DATA ---
typedef struct {
    int rod_length;
    int *prices;    // Array of prices for lengths 1 to rod_length
    int *revenue;   // Array for DP results (tabulation)
    int final_result; // Final computed result
} BenchmarkData;

static BenchmarkData g_data;

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <rod_length> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.rod_length = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if (g_data.rod_length <= 0) {
        fprintf(stderr, "FATAL: rod_length must be a positive integer.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate memory. Prices are for lengths 1...rod_length.
    // Index 0 is unused for prices but makes indexing easier.
    g_data.prices = (int *)malloc(sizeof(int) * (g_data.rod_length + 1));
    if (g_data.prices == NULL) {
        fprintf(stderr, "FATAL: Memory allocation for prices failed.\n");
        exit(1);
    }

    // Allocate memory for the DP table.
    g_data.revenue = (int *)malloc(sizeof(int) * (g_data.rod_length + 1));
    if (g_data.revenue == NULL) {
        fprintf(stderr, "FATAL: Memory allocation for revenue table failed.\n");
        free(g_data.prices);
        exit(1);
    }
    
    // Generate random prices for each possible cut length.
    // Price for a piece of length i is stored at prices[i].
    g_data.prices[0] = 0; // Price for length 0 is 0.
    for (int i = 1; i <= g_data.rod_length; i++) {
        // Generate a plausible price, capped to keep sums reasonable.
        g_data.prices[i] = (mt_rand() % (i * 5 + 50)) + 1;
    }
}

void run_computation() {
    int n = g_data.rod_length;
    int *revenue = g_data.revenue;
    int *prices = g_data.prices;

    revenue[0] = 0;

    // Build the revenue table from the bottom up using tabulation.
    // For each rod length i from 1 to n...
    for (int i = 1; i <= n; i++) {
        int max_val = -1;
        // ...consider all possible first cuts of length j (where j is from 1 to i).
        for (int j = 1; j <= i; j++) {
            // The value is the price of the first cut piece (length j) plus
            // the max revenue from the remaining piece (length i-j).
            int current_val = prices[j] + revenue[i - j];
            if (current_val > max_val) {
                max_val = current_val;
            }
        }
        revenue[i] = max_val;
    }

    g_data.final_result = revenue[n];
}

void cleanup() {
    free(g_data.prices);
    free(g_data.revenue);
    g_data.prices = NULL;
    g_data.revenue = NULL;
}

// --- MAIN FUNCTION ---
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final result to stdout
    printf("%d\n", g_data.final_result);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
