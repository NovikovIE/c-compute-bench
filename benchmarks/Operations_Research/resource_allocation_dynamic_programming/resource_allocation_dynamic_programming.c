/**
 * @file resource_allocation_dynamic_programming.c
 * @brief A benchmark simulating a resource allocation problem using dynamic programming.
 * 
 * This program models the classic 0/1 knapsack problem, a common task in operations
 * research. Given a set of projects, each with a cost and a potential return, and a
 * total budget, the goal is to select a subset of projects that maximizes the total
 * return without exceeding the budget.
 *
 * The problem is solved using a space-optimized dynamic programming approach.
 * The time complexity is O(num_projects * total_budget), and the space complexity is
 * O(total_budget).
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- Mersenne Twister (MT19937) PRNG --- (DO NOT MODIFY)
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
// --- End of MT19937 --- C

// --- Benchmark Globals ---

// Parameters
int num_projects;
int total_budget;

// Data structures
int *project_costs;
int *project_returns;

// DP table for computation
int *dp_table;

// Final result
int final_result = 0;

// --- Utility Functions ---
#ifndef max
    #define max(a,b) (((a) > (b)) ? (a) : (b))
#endif

// --- Benchmark Functions ---

/**
 * @brief Parses arguments, allocates memory, and generates random project data.
 */
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_projects> <total_budget> <seed>\n", argv[0]);
        exit(1);
    }

    num_projects = atoi(argv[1]);
    total_budget = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (num_projects <= 0 || total_budget <= 0) {
        fprintf(stderr, "FATAL: num_projects and total_budget must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate memory for project data
    project_costs = (int*)malloc(num_projects * sizeof(int));
    project_returns = (int*)malloc(num_projects * sizeof(int));
    if (!project_costs || !project_returns) {
        fprintf(stderr, "FATAL: Memory allocation failed for project data.\n");
        exit(1);
    }

    // Generate random data for each project
    for (int i = 0; i < num_projects; i++) {
        // Costs are relatively low to allow for more combinations
        project_costs[i] = 1 + (mt_rand() % 200);
        // Returns are higher to make the sum significant
        project_returns[i] = 1 + (mt_rand() % 5000);
    }

    // Allocate and initialize the DP table. calloc initializes memory to zero.
    dp_table = (int*)calloc(total_budget + 1, sizeof(int));
    if (!dp_table) {
        fprintf(stderr, "FATAL: Memory allocation failed for DP table.\n");
        exit(1);
    }
}

/**
 * @brief Runs the core dynamic programming algorithm.
 * 
 * Solves the 0/1 knapsack problem using a space-optimized DP table.
 * For each project, it iterates through the budget from max to the project's cost,
 * updating the maximum possible return for each budget point.
 */
void run_computation() {
    for (int i = 0; i < num_projects; i++) {
        int cost = project_costs[i];
        int ret = project_returns[i];
        for (int b = total_budget; b >= cost; b--) {
            dp_table[b] = max(dp_table[b], ret + dp_table[b - cost]);
        }
    }
    final_result = dp_table[total_budget];
}

/**
 * @brief Frees all memory allocated in setup_benchmark.
 */
void cleanup() {
    free(project_costs);
    free(project_returns);
    free(dp_table);
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

    // Print final result to stdout
    printf("%d\n", final_result);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
