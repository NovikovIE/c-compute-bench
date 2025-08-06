#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

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

// A constraint is a pair of variables that must have different values.
typedef struct {
    int var1;
    int var2;
} Constraint;

// Global struct to hold all benchmark data
struct {
    int num_variables;
    int domain_size_per_variable;
    int max_steps;
    int num_constraints;

    int *assignment; // current value for each variable
    Constraint *constraints; // list of all constraints

    // Adjacency list representation of the constraint graph
    int **adj_list; 
    int *adj_list_sizes;

    // Workspace arrays for computation
    int *conflicted_vars; 
    int *conflict_count; 

    int final_result; // The result to be printed
} *g_data;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <num_vars> <domain_size> <max_steps> <num_constraints> <seed>\n", argv[0]);
        exit(1);
    }

    g_data = (void *)malloc(sizeof(*g_data));
    if (!g_data) { perror("malloc failed"); exit(1); }

    g_data->num_variables = atoi(argv[1]);
    g_data->domain_size_per_variable = atoi(argv[2]);
    g_data->max_steps = atoi(argv[3]);
    g_data->num_constraints = atoi(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    mt_seed(seed);

    // Allocate memory
    g_data->assignment = (int *)malloc(g_data->num_variables * sizeof(int));
    g_data->constraints = (Constraint *)malloc(g_data->num_constraints * sizeof(Constraint));
    g_data->adj_list_sizes = (int *)calloc(g_data->num_variables, sizeof(int));
    g_data->adj_list = (int **)malloc(g_data->num_variables * sizeof(int *));
    g_data->conflicted_vars = (int *)malloc(g_data->num_variables * sizeof(int));
    g_data->conflict_count = (int *)calloc(g_data->num_variables, sizeof(int));

    if (!g_data->assignment || !g_data->constraints || !g_data->adj_list_sizes || 
        !g_data->adj_list || !g_data->conflicted_vars || !g_data->conflict_count) {
        perror("malloc failed for data arrays");
        exit(1);
    }

    // Generate random constraints (var1 != var2, binary inequality)
    for (int i = 0; i < g_data->num_constraints; ++i) {
        int v1 = mt_rand() % g_data->num_variables;
        int v2;
        do {
            v2 = mt_rand() % g_data->num_variables;
        } while (v1 == v2);
        g_data->constraints[i] = (Constraint){v1, v2};
        g_data->adj_list_sizes[v1]++;
        g_data->adj_list_sizes[v2]++;
    }

    // Build adjacency list from constraints
    for (int i = 0; i < g_data->num_variables; ++i) {
        g_data->adj_list[i] = (int *)malloc(g_data->adj_list_sizes[i] * sizeof(int));
        if (!g_data->adj_list[i]) { perror("malloc failed for adj_list row"); exit(1); }
        g_data->adj_list_sizes[i] = 0; // reset to use as a fill counter
    }
    for (int i = 0; i < g_data->num_constraints; ++i) {
        int v1 = g_data->constraints[i].var1;
        int v2 = g_data->constraints[i].var2;
        g_data->adj_list[v1][g_data->adj_list_sizes[v1]++] = v2;
        g_data->adj_list[v2][g_data->adj_list_sizes[v2]++] = v1;
    }

    // Generate initial random assignment
    for (int i = 0; i < g_data->num_variables; ++i) {
        g_data->assignment[i] = mt_rand() % g_data->domain_size_per_variable;
    }

    // Calculate initial conflict counts
    for (int i = 0; i < g_data->num_constraints; ++i) {
        int v1 = g_data->constraints[i].var1;
        int v2 = g_data->constraints[i].var2;
        if (g_data->assignment[v1] == g_data->assignment[v2]) {
            g_data->conflict_count[v1]++;
            g_data->conflict_count[v2]++;
        }
    }
}

void run_computation() {
    int total_reassignments = 0;
    int* best_values = (int*)malloc(g_data->domain_size_per_variable * sizeof(int));

    for (int step = 0; step < g_data->max_steps; ++step) {
        int num_conflicted = 0;
        for (int i = 0; i < g_data->num_variables; ++i) {
            if (g_data->conflict_count[i] > 0) {
                g_data->conflicted_vars[num_conflicted++] = i;
            }
        }

        if (num_conflicted == 0) {
            break; // Solution found
        }

        int var_to_change = g_data->conflicted_vars[mt_rand() % num_conflicted];
        int min_conflicts = g_data->num_variables + 1; // more than max possible conflicts
        int num_best = 0;

        // Find the value(s) that minimize conflicts for var_to_change
        for (int val = 0; val < g_data->domain_size_per_variable; ++val) {
            int current_conflicts = 0;
            for (int i = 0; i < g_data->adj_list_sizes[var_to_change]; ++i) {
                int neighbor = g_data->adj_list[var_to_change][i];
                if (g_data->assignment[neighbor] == val) {
                    current_conflicts++;
                }
            }
            if (current_conflicts < min_conflicts) {
                min_conflicts = current_conflicts;
                best_values[0] = val;
                num_best = 1;
            } else if (current_conflicts == min_conflicts) {
                best_values[num_best++] = val;
            }
        }
        
        int best_value = best_values[mt_rand() % num_best];
        int old_value = g_data->assignment[var_to_change];

        if (old_value != best_value) {
            total_reassignments++;

            // Update conflict counts for neighbors
            for (int i = 0; i < g_data->adj_list_sizes[var_to_change]; ++i) {
                int neighbor = g_data->adj_list[var_to_change][i];
                if (g_data->assignment[neighbor] == old_value) {
                    g_data->conflict_count[neighbor]--;
                }
                if (g_data->assignment[neighbor] == best_value) {
                    g_data->conflict_count[neighbor]++;
                }
            }
            g_data->assignment[var_to_change] = best_value;
            g_data->conflict_count[var_to_change] = min_conflicts;
        }
    }

    g_data->final_result = total_reassignments;
    free(best_values);
}

void cleanup() {
    for (int i = 0; i < g_data->num_variables; ++i) {
        free(g_data->adj_list[i]);
    }
    free(g_data->adj_list);
    free(g_data->adj_list_sizes);
    free(g_data->constraints);
    free(g_data->assignment);
    free(g_data->conflicted_vars);
    free(g_data->conflict_count);
    g_data->final_result = 0;
    free(g_data);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
    printf("%d\n", g_data->final_result);
    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}