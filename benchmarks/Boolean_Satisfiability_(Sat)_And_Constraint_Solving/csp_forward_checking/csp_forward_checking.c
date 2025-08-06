#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>

// --- BEGIN: Mersenne Twister (DO NOT MODIFY) ---
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
// --- END: Mersenne Twister ---

// --- Benchmark Globals ---
typedef struct {
    int num_variables;
    int domain_size;
    int tightness_percentage;

    // Constraint representation: true = DISALLOWED pair
    bool *constraints;
    int **constraint_indices; // Maps (i, j) where i < j to a constraint table index

    // Solver state
    int *assignment;
    bool **current_domains; // Domains pruned during search

    long long solution_count; // Final result
} CspProblem;

CspProblem G;

// --- Function Prototypes ---
void setup_benchmark(int argc, char *argv[]);
void run_computation();
void cleanup();
long long solve_recursive(int var_index);

// --- Main and Timing ---
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%lld\n", G.solution_count);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}

// --- Benchmark Implementation ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_variables domain_size_per_variable constraint_tightness_percentage seed\n", argv[0]);
        exit(1);
    }

    G.num_variables = atoi(argv[1]);
    G.domain_size = atoi(argv[2]);
    G.tightness_percentage = atoi(argv[3]);
    mt_seed(atoi(argv[4]));

    // Allocate memory
    G.assignment = (int *)malloc(G.num_variables * sizeof(int));
    G.current_domains = (bool **)malloc(G.num_variables * sizeof(bool *));
    for (int i = 0; i < G.num_variables; i++) {
        G.current_domains[i] = (bool *)malloc(G.domain_size * sizeof(bool));
        for (int j = 0; j < G.domain_size; j++) {
            G.current_domains[i][j] = true; // All values initially available
        }
    }

    int num_constraints = G.num_variables * (G.num_variables - 1) / 2;
    int constraint_table_size = G.domain_size * G.domain_size;
    G.constraints = (bool *)malloc((size_t)num_constraints * constraint_table_size * sizeof(bool));
    memset(G.constraints, 0, (size_t)num_constraints * constraint_table_size * sizeof(bool)); // All pairs allowed by default

    G.constraint_indices = (int **)malloc(G.num_variables * sizeof(int *));
    for (int i = 0; i < G.num_variables; i++) {
        G.constraint_indices[i] = (int *)malloc(G.num_variables * sizeof(int));
    }

    // Generate constraints
    int num_disallowed_pairs = (int)(G.tightness_percentage / 100.0 * constraint_table_size);
    int* pair_indices_shuffler = (int*)malloc(constraint_table_size * sizeof(int));
    for(int i = 0; i < constraint_table_size; ++i) pair_indices_shuffler[i] = i;

    int constraint_idx_counter = 0;
    for (int i = 0; i < G.num_variables; i++) {
        for (int j = i + 1; j < G.num_variables; j++) {
            G.constraint_indices[i][j] = constraint_idx_counter;
            long long constraint_base_offset = (long long)constraint_idx_counter * constraint_table_size;

            // Shuffle pair indices to select random pairs to disallow
            for (int k = constraint_table_size - 1; k > 0; k--) {
                int l = mt_rand() % (k + 1);
                int temp = pair_indices_shuffler[k];
                pair_indices_shuffler[k] = pair_indices_shuffler[l];
                pair_indices_shuffler[l] = temp;
            }
            
            // Mark pairs as disallowed
            for (int k = 0; k < num_disallowed_pairs; k++) {
                G.constraints[constraint_base_offset + pair_indices_shuffler[k]] = true;
            }            
            constraint_idx_counter++;
        }
    }
    free(pair_indices_shuffler);
}


long long solve_recursive(int var_index) {
    if (var_index == G.num_variables) {
        return 1;
    }

    long long count = 0;

    for (int v = 0; v < G.domain_size; ++v) {
        if (G.current_domains[var_index][v]) {
            G.assignment[var_index] = v;

            // Store prunings to be undone upon backtracking
            int max_pruned = (G.num_variables - (var_index + 1)) * G.domain_size;
            int *pruned_vars = (int *)malloc(max_pruned * sizeof(int));
            int *pruned_vals = (int *)malloc(max_pruned * sizeof(int));
            int num_pruned = 0;
            bool consistent_path = true;

            // --- Forward Check ---
            for (int future_var = var_index + 1; future_var < G.num_variables; ++future_var) {
                long long c_base_idx = (long long)G.constraint_indices[var_index][future_var] * G.domain_size * G.domain_size;
                for (int future_val = 0; future_val < G.domain_size; ++future_val) {
                    if (G.current_domains[future_var][future_val]) {
                        long long flat_idx = c_base_idx + (long long)v * G.domain_size + future_val;
                        if (G.constraints[flat_idx]) { // if pair is disallowed
                            G.current_domains[future_var][future_val] = false;
                            pruned_vars[num_pruned] = future_var;
                            pruned_vals[num_pruned] = future_val;
                            num_pruned++;
                        }
                    }
                }

                // Check if domain of future_var became empty
                bool domain_wiped_out = true;
                for (int i = 0; i < G.domain_size; ++i) {
                    if (G.current_domains[future_var][i]) {
                        domain_wiped_out = false;
                        break;
                    }
                }
                if (domain_wiped_out) {
                    consistent_path = false;
                    break;
                }
            }

            if (consistent_path) {
                count += solve_recursive(var_index + 1);
            }

            // --- Backtrack: Undo Pruning ---
            for (int i = 0; i < num_pruned; ++i) {
                G.current_domains[pruned_vars[i]][pruned_vals[i]] = true;
            }
            free(pruned_vars);
            free(pruned_vals);
        }
    }

    return count;
}

void run_computation() {
    G.solution_count = solve_recursive(0);
}

void cleanup() {
    free(G.assignment);
    free(G.constraints);
    
    for (int i = 0; i < G.num_variables; i++) {
        free(G.current_domains[i]);
        free(G.constraint_indices[i]);
    }
    free(G.current_domains);
    free(G.constraint_indices);
}
