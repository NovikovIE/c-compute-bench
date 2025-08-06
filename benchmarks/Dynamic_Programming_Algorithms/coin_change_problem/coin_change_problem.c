#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <limits.h>

// --- BEGIN MERSENNE TWISTER (MT19937) --- Do Not Modify ---
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

// Global data structure to hold benchmark data
typedef struct {
    int amount;
    int num_coin_types;
    int *coins;
    int *dp_table;
    int final_result;
} BenchmarkData;

BenchmarkData g_data;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <amount> <num_coin_types> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.amount = atoi(argv[1]);
    g_data.num_coin_types = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    mt_seed(seed);

    // Allocate memory for coin denominations and the DP table
    g_data.coins = (int *)malloc(g_data.num_coin_types * sizeof(int));
    if (g_data.coins == NULL) {
        fprintf(stderr, "Failed to allocate memory for coins.\n");
        exit(1);
    }

    g_data.dp_table = (int *)malloc((g_data.amount + 1) * sizeof(int));
    if (g_data.dp_table == NULL) {
        fprintf(stderr, "Failed to allocate memory for DP table.\n");
        free(g_data.coins);
        exit(1);
    }

    // Generate unique coin denominations
    // We add coin '1' to ensure a solution always exists.
    g_data.coins[0] = 1;
    int coin_max_val = g_data.amount > 1000 ? g_data.amount / 10 : 100;
    for (int i = 1; i < g_data.num_coin_types; ++i) {
        int new_coin;
        int is_duplicate;
        do {
            is_duplicate = 0;
            new_coin = 2 + (mt_rand() % coin_max_val);
            for (int j = 0; j < i; ++j) {
                if (g_data.coins[j] == new_coin) {
                    is_duplicate = 1;
                    break;
                }
            }
        } while (is_duplicate);
        g_data.coins[i] = new_coin;
    }

    // Initialize the DP table. amount + 1 serves as infinity.
    int infinity = g_data.amount + 1;
    g_data.dp_table[0] = 0;
    for (int i = 1; i <= g_data.amount; ++i) {
        g_data.dp_table[i] = infinity;
    }
    
    g_data.final_result = 0;
}

void run_computation() {
    int infinity = g_data.amount + 1;

    // Bottom-up DP (Tabulation)
    for (int i = 0; i < g_data.num_coin_types; i++) {
        int coin = g_data.coins[i];
        for (int j = coin; j <= g_data.amount; j++) {
            if (g_data.dp_table[j - coin] != infinity) {
                int current_ways = g_data.dp_table[j];
                int new_ways = g_data.dp_table[j - coin] + 1;
                if (new_ways < current_ways) {
                    g_data.dp_table[j] = new_ways;
                }
            }
        }
    }

    g_data.final_result = (g_data.dp_table[g_data.amount] == infinity) ? -1 : g_data.dp_table[g_data.amount];
}

void cleanup() {
    free(g_data.coins);
    free(g_data.dp_table);
}

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

    // Print the timing info to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
