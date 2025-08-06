#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <limits.h>

// --- BEGIN MERSENNE TWISTER (DO NOT MODIFY) ---
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
// --- END MERSENNE TWISTER ---

// --- BENCHMARK DATA AND PARAMETERS ---
typedef struct {
    int num_time_periods;
    int max_inventory_level;
    long long *demand;
    long long **dp_table;
    long long final_result;

    // Cost parameters for the inventory model
    const long long production_cost;
    const long long holding_cost;
    const long long fixed_ordering_cost;
} BenchmarkData;

BenchmarkData g_data = {0, 0, NULL, NULL, 0, 5, 1, 100};

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_time_periods> <max_inventory_level> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_time_periods = atoi(argv[1]);
    g_data.max_inventory_level = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    if (g_data.num_time_periods <= 0 || g_data.max_inventory_level <= 0) {
        fprintf(stderr, "Parameters must be positive integers.\n");
        exit(1);
    }

    mt_seed(seed);

    g_data.demand = (long long*)malloc(g_data.num_time_periods * sizeof(long long));
    if (g_data.demand == NULL) {
        fprintf(stderr, "Failed to allocate memory for demand array.\n");
        exit(1);
    }

    int M = g_data.max_inventory_level;
    g_data.dp_table = (long long**)malloc(g_data.num_time_periods * sizeof(long long*));
    if (g_data.dp_table == NULL) {
        fprintf(stderr, "Failed to allocate memory for dp_table.\n");
        free(g_data.demand);
        exit(1);
    }
    for (int i = 0; i < g_data.num_time_periods; i++) {
        g_data.dp_table[i] = (long long*)malloc((M + 1) * sizeof(long long));
        if (g_data.dp_table[i] == NULL) {
            fprintf(stderr, "Failed to allocate memory for dp_table row.\n");
            // Cleanup already allocated memory
            for (int k = 0; k < i; k++) free(g_data.dp_table[k]);
            free(g_data.dp_table);
            free(g_data.demand);
            exit(1);
        }
    }

    // Generate random demand for each period, ensuring it's not excessively large.
    for (int t = 0; t < g_data.num_time_periods; t++) {
        g_data.demand[t] = (mt_rand() % (M / 4)) + (M / 10);
    }
}

void run_computation() {
    int N = g_data.num_time_periods;
    int M = g_data.max_inventory_level;

    // Base case: t = 0
    for (int i = 0; i <= M; i++) { // i is ending inventory for period 0
        long long production = g_data.demand[0] + i;
        long long cost = g_data.production_cost * production + g_data.holding_cost * i;
        if (production > 0) {
            cost += g_data.fixed_ordering_cost;
        }
        g_data.dp_table[0][i] = cost;
    }

    // Dynamic Programming recurrence: t = 1 to N-1
    for (int t = 1; t < N; t++) {
        for (int j = 0; j <= M; j++) { // j is ending inventory for period t
            g_data.dp_table[t][j] = LLONG_MAX;
            for (int i = 0; i <= M; i++) { // i is ending inventory for period t-1

                if (g_data.dp_table[t - 1][i] == LLONG_MAX) continue;

                long long production = (long long)j - i + g_data.demand[t];
                if (production < 0) continue; // This transition is not feasible

                long long period_cost = g_data.production_cost * production + g_data.holding_cost * j;
                if (production > 0) {
                    period_cost += g_data.fixed_ordering_cost;
                }

                long long new_total_cost = g_data.dp_table[t - 1][i] + period_cost;
                
                if (new_total_cost < g_data.dp_table[t][j]) {
                    g_data.dp_table[t][j] = new_total_cost;
                }
            }
        }
    }

    // Find the minimum total cost at the end of the planning horizon
    long long min_cost = LLONG_MAX;
    for (int i = 0; i <= M; i++) {
        if (g_data.dp_table[N - 1][i] < min_cost) {
            min_cost = g_data.dp_table[N - 1][i];
        }
    }
    g_data.final_result = min_cost;
}

void cleanup() {
    if (g_data.dp_table != NULL) {
        for (int i = 0; i < g_data.num_time_periods; i++) {
            free(g_data.dp_table[i]);
        }
        free(g_data.dp_table);
    }
    if (g_data.demand != NULL) {
        free(g_data.demand);
    }
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
    printf("%lld\n", g_data.final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
