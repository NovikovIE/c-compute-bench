#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <limits.h>

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
// END of Mersenne Twister

// Benchmark-specific global variables
static int num_eggs;
static int num_floors;
static int** dp_table;
static int final_result;

// Helper to find maximum of two integers
static inline int max(int a, int b) {
    return a > b ? a : b;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_eggs> <num_floors> <seed>\n", argv[0]);
        exit(1);
    }

    num_eggs = atoi(argv[1]);
    num_floors = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (num_eggs <= 0 || num_floors < 0) {
        fprintf(stderr, "Error: num_eggs must be > 0 and num_floors must be >= 0.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate DP table: (num_eggs + 1) x (num_floors + 1)
    dp_table = (int **)malloc((num_eggs + 1) * sizeof(int *));
    if (!dp_table) {
        fprintf(stderr, "Memory allocation failed for dp_table rows\n");
        exit(1);
    }
    for (int i = 0; i <= num_eggs; i++) {
        dp_table[i] = (int *)malloc((num_floors + 1) * sizeof(int));
        if (!dp_table[i]) {
            fprintf(stderr, "Memory allocation failed for dp_table columns\n");
            for(int j = 0; j < i; j++) free(dp_table[j]);
            free(dp_table);
            exit(1);
        }
    }
}

void run_computation() {
    // Base case: 1 egg, f floors -> f trials
    for (int f = 1; f <= num_floors; f++) {
        dp_table[1][f] = f;
    }

    // Base cases: e eggs, 0 floors -> 0 trials; e eggs, 1 floor -> 1 trial
    for (int e = 1; e <= num_eggs; e++) {
        dp_table[e][0] = 0;
        dp_table[e][1] = 1;
    }

    // Fill the rest of the table using DP
    for (int e = 2; e <= num_eggs; e++) {
        for (int f = 2; f <= num_floors; f++) {
            dp_table[e][f] = INT_MAX;
            for (int x = 1; x <= f; x++) {
                // If we drop from floor x:
                // 1. Egg breaks: we have e-1 eggs, and need to check floors 1 to x-1.
                // 2. Egg survives: we have e eggs, and need to check floors x+1 to f (f-x floors).
                int res = 1 + max(dp_table[e - 1][x - 1], dp_table[e][f - x]);
                if (res < dp_table[e][f]) {
                    dp_table[e][f] = res;
                }
            }
        }
    }

    final_result = dp_table[num_eggs][num_floors];
}

void cleanup() {
    for (int i = 0; i <= num_eggs; i++) {
        free(dp_table[i]);
    }
    free(dp_table);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%d\n", final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
