#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

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

// Parameters
int g_num_variables;
int g_domain_size;
int g_tightness_pct;

// Data structures for the Constraint Satisfaction Problem (CSP)
// Using char as a boolean (0 or 1) to save space.
// Represents a dense table of binary constraints.
char* g_constraints;
// Stores the current assignment of values to variables during search.
int* g_assignment;
// Final result: the total number of valid solutions found.
long long g_solution_count;

// --- HELPER FUNCTIONS ---

// Calculates the index into the flattened 4D constraints array.
// Constraint(var1, var2, val1, val2)
static inline size_t get_constraint_idx(int var1, int var2, int val1, int val2) {
    size_t idx = 0;
    idx += (size_t)var1 * g_num_variables * g_domain_size * g_domain_size;
    idx += (size_t)var2 * g_domain_size * g_domain_size;
    idx += (size_t)val1 * g_domain_size;
    idx += val2;
    return idx;
}

// Checks if the current assignment for `var_idx` is consistent with all
// previously assigned variables.
static int is_consistent(int var_idx) {
    int current_value = g_assignment[var_idx];
    for (int i = 0; i < var_idx; i++) {
        int previous_value = g_assignment[i];
        // Look up the constraint from our table
        size_t idx = get_constraint_idx(var_idx, i, current_value, previous_value);
        if (!g_constraints[idx]) {
            return 0; // The assignment violates a constraint
        }
    }
    return 1; // Consistent with all previous assignments
}

// Recursive backtracking function to find and count all solutions.
static long long backtrack_search(int var_idx) {
    // Base case: If all variables have been assigned, we have found a solution.
    if (var_idx == g_num_variables) {
        return 1;
    }

    long long count = 0;
    // Try assigning each possible value from the domain to the current variable.
    for (int val = 0; val < g_domain_size; val++) {
        g_assignment[var_idx] = val;
        // If the assignment is consistent with previous ones, recurse.
        if (is_consistent(var_idx)) {
            count += backtrack_search(var_idx + 1);
        }
    }
    return count;
}

// --- BENCHMARK CORE FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_variables> <domain_size> <tightness_pct> <seed>\n", argv[0]);
        exit(1);
    }

    g_num_variables = atoi(argv[1]);
    g_domain_size = atoi(argv[2]);
    g_tightness_pct = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    mt_seed(seed);

    if (g_num_variables <= 0 || g_domain_size <= 0 || g_tightness_pct < 0 || g_tightness_pct > 100) {
        fprintf(stderr, "Invalid parameters.\n");
        exit(1);
    }

    // Allocate memory for a dense N x N x D x D constraint table.
    size_t constraint_table_size = (size_t)g_num_variables * g_num_variables * g_domain_size * g_domain_size;
    g_constraints = (char*)malloc(constraint_table_size * sizeof(char));
    if (!g_constraints) {
        fprintf(stderr, "Failed to allocate memory for constraints.\n");
        exit(1);
    }

    // Generate random constraints.
    for (int i = 0; i < g_num_variables; i++) {
        for (int j = i; j < g_num_variables; j++) {
            // No constraints for a variable with itself.
            if (i == j) continue;

            for (int val_i = 0; val_i < g_domain_size; val_i++) {
                for (int val_j = 0; val_j < g_domain_size; val_j++) {
                    // A pair of values is allowed if a random number is outside the tightness percentage.
                    char is_allowed = (mt_rand() % 100) >= g_tightness_pct;
                    g_constraints[get_constraint_idx(i, j, val_i, val_j)] = is_allowed;
                    // Make constraints symmetric for easier lookup.
                    g_constraints[get_constraint_idx(j, i, val_j, val_i)] = is_allowed;
                }
            }
        }
    }

    // Allocate memory for the assignment array.
    g_assignment = (int*)malloc(g_num_variables * sizeof(int));
    if (!g_assignment) {
        fprintf(stderr, "Failed to allocate memory for assignment array.\n");
        free(g_constraints);
        exit(1);
    }
}

void run_computation() {
    g_solution_count = backtrack_search(0);
}

void cleanup() {
    free(g_constraints);
    free(g_assignment);
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

    // Print final result (total number of solutions) to stdout.
    printf("%lld\n", g_solution_count);

    // Print timing information to stderr.
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}