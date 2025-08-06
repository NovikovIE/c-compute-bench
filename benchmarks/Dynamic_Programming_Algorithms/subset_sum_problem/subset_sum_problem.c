#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- Mersenne Twister (MT19937) --- (DO NOT MODIFY)
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
typedef struct {
    int set_size;
    int target_sum;
    uint32_t *set;
    char **dp_table; // Using char as boolean (0 or 1)
    int result;
} BenchmarkData;

BenchmarkData g_data;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "FATAL: Usage: %s <set_size> <target_sum> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.set_size = atoi(argv[1]);
    g_data.target_sum = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    if (g_data.set_size <= 0 || g_data.target_sum < 0) {
        fprintf(stderr, "FATAL: set_size must be positive and target_sum must be non-negative.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate set and fill with random numbers
    g_data.set = (uint32_t *)malloc(g_data.set_size * sizeof(uint32_t));
    if (g_data.set == NULL) {
        fprintf(stderr, "FATAL: Failed to allocate memory for the set.\n");
        exit(1);
    }
    for (int i = 0; i < g_data.set_size; i++) {
        g_data.set[i] = (mt_rand() % 500) + 1; // Numbers from 1 to 500
    }

    // Allocate the DP table
    g_data.dp_table = (char **)malloc((g_data.set_size + 1) * sizeof(char *));
    if (g_data.dp_table == NULL) {
        fprintf(stderr, "FATAL: Failed to allocate memory for DP table rows.\n");
        free(g_data.set);
        exit(1);
    }

    for (int i = 0; i <= g_data.set_size; i++) {
        g_data.dp_table[i] = (char *)malloc((g_data.target_sum + 1) * sizeof(char));
        if (g_data.dp_table[i] == NULL) {
            fprintf(stderr, "FATAL: Failed to allocate memory for DP table columns.\n");
            // Cleanup already allocated memory
            for (int k = 0; k < i; k++) {
                free(g_data.dp_table[k]);
            }
            free(g_data.dp_table);
            free(g_data.set);
            exit(1);
        }
    }
    g_data.result = 0;
}

void run_computation() {
    // dp_table[i][j] will be true if sum j can be obtained using a subset of first i items

    // Base Case 1: If sum is 0, the answer is always true (empty subset)
    for (int i = 0; i <= g_data.set_size; i++) {
        g_data.dp_table[i][0] = 1;
    }

    // Base Case 2: If the set is empty (i=0), no sum > 0 can be obtained
    for (int j = 1; j <= g_data.target_sum; j++) {
        g_data.dp_table[0][j] = 0;
    }

    // Fill the rest of the DP table using tabulation
    for (int i = 1; i <= g_data.set_size; i++) {
        for (int j = 1; j <= g_data.target_sum; j++) {
            // If the i-th element is not included
            g_data.dp_table[i][j] = g_data.dp_table[i - 1][j];

            // If the i-th element can be included
            uint32_t current_item = g_data.set[i - 1];
            if (j >= current_item) {
                g_data.dp_table[i][j] = g_data.dp_table[i][j] || g_data.dp_table[i - 1][j - current_item];
            }
        }
    }

    g_data.result = g_data.dp_table[g_data.set_size][g_data.target_sum];
}

void cleanup() {
    if (g_data.dp_table) {
        for (int i = 0; i <= g_data.set_size; i++) {
            free(g_data.dp_table[i]);
        }
        free(g_data.dp_table);
    }
    free(g_data.set);
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

    // Print final result to stdout
    printf("%d\n", g_data.result);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
