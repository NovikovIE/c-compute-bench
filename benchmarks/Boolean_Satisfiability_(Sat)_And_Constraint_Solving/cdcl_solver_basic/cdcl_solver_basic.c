#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

// --- Mersenne Twister (MT19937) Generator - Do Not Modify ---
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
// --- End of Mersenne Twister ---

// A literal is a variable (positive int) or its negation (negative int).
// E.g., variable x3 is represented as 3, and its negation !x3 as -3.
// Variables are 1-indexed.

typedef struct {
    int* literals;
    int size;
} Clause;

typedef struct {
    // Problem specification
    int num_variables;
    int num_clauses;
    int literals_per_clause; // From clause_to_variable_ratio
    Clause* clauses;

    // Solver state
    int* assignment;  // 0: unassigned, 1: TRUE, -1: FALSE (index i for var i)
    int* decision_level; // Decision level for each variable assignment
    int* trail;       // Order of variable assignments
    int trail_count;
    int decision_level_current;
    
    // Result accumulator
    long long conflicts_found;
} BenchmarkData;

BenchmarkData* g_data = NULL;

// --- Function Prototypes ---
void setup_benchmark(int argc, char *argv[]);
void run_computation();
void cleanup();

// --- Main Function ---
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%lld\n", g_data->conflicts_found);
    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_variables num_clauses clause_to_variable_ratio seed\n", argv[0]);
        exit(1);
    }

    g_data = (BenchmarkData*)malloc(sizeof(BenchmarkData));
    if (!g_data) {
        perror("Failed to allocate memory for benchmark data");
        exit(1);
    }

    g_data->num_variables = atoi(argv[1]);
    g_data->num_clauses = atoi(argv[2]);
    // 'clause_to_variable_ratio' is interpreted as 'k' in k-SAT.
    g_data->literals_per_clause = (int)round(atof(argv[3]));
    uint32_t seed = (uint32_t)atoi(argv[4]);
    mt_seed(seed);

    int k = g_data->literals_per_clause;
    if (k < 1) {
        fprintf(stderr, "Literals per clause must be at least 1.\n");
        exit(1);
    }

    // Generate a random k-SAT formula
    g_data->clauses = (Clause*)malloc(g_data->num_clauses * sizeof(Clause));
    bool* var_is_used = (bool*)malloc((g_data->num_variables + 1) * sizeof(bool)); 

    for (int i = 0; i < g_data->num_clauses; ++i) {
        g_data->clauses[i].size = k;
        g_data->clauses[i].literals = (int*)malloc(k * sizeof(int));
        memset(var_is_used, 0, (g_data->num_variables + 1) * sizeof(bool));

        for (int j = 0; j < k; ) {
            int var = (mt_rand() % g_data->num_variables) + 1;
            if (!var_is_used[var]) {
                var_is_used[var] = true;
                int sign = (mt_rand() % 2) == 0 ? 1 : -1;
                g_data->clauses[i].literals[j] = var * sign;
                j++;
            }
        }
    }
    free(var_is_used);

    // Initialize solver state
    // Use (num_variables + 1) for 1-based indexing
    g_data->assignment = (int*)calloc(g_data->num_variables + 1, sizeof(int));
    g_data->decision_level = (int*)calloc(g_data->num_variables + 1, sizeof(int));
    g_data->trail = (int*)malloc((g_data->num_variables + 1) * sizeof(int));
    g_data->trail_count = 0;
    g_data->decision_level_current = 0;
    g_data->conflicts_found = 0;
}

/**
 * Implements a basic iterative DPLL (Davis-Putnam-Logemann-Loveland) algorithm.
 * This is a foundational algorithm for modern SAT solvers.
 * The main loop consists of propagation, decision, and backtracking.
 */
void run_computation() {
    int* assignment = g_data->assignment;
    int* decision_level = g_data->decision_level;
    int* trail = g_data->trail;

    int propagation_head = 0;

    while (g_data->trail_count < g_data->num_variables) {
        // 1. PROPAGATE: Check for unit clauses resulting from recent assignments.
        // A full solver uses watch lists for efficiency; here we scan all clauses.
        if (propagation_head < g_data->trail_count) {
            for (int i = 0; i < g_data->num_clauses; ++i) {
                Clause* c = &g_data->clauses[i];
                int unassigned_count = 0;
                int unassigned_lit = 0;
                bool is_satisfied = false;

                for (int j = 0; j < c->size; ++j) {
                    int lit = c->literals[j];
                    int var = abs(lit);
                    int val = lit > 0 ? 1 : -1;

                    if (assignment[var] == val) {
                        is_satisfied = true;
                        break;
                    }
                    if (assignment[var] == 0) {
                        unassigned_count++;
                        unassigned_lit = lit;
                    }
                }

                if (is_satisfied) continue;

                if (unassigned_count == 0) { // CONFLICT
                    goto handle_conflict;
                }

                if (unassigned_count == 1) { // UNIT CLAUSE
                    int var = abs(unassigned_lit);
                    if (assignment[var] == 0) {
                        assignment[var] = unassigned_lit > 0 ? 1 : -1;
                        decision_level[var] = g_data->decision_level_current;
                        trail[g_data->trail_count++] = var;
                    }
                }
            }
            propagation_head = g_data->trail_count;
            continue; // Re-evaluate loop condition after potential propagation
        }

        // 2. DECIDE: If no conflict, pick an unassigned variable and assign it.
        g_data->decision_level_current++;
        int next_var = 0;
        for (int i = 1; i <= g_data->num_variables; i++) {
            if (assignment[i] == 0) {
                next_var = i;
                break;
            }
        }

        if (next_var == 0) break; // SATISFIABLE

        assignment[next_var] = 1; // Decide TRUE
        decision_level[next_var] = g_data->decision_level_current;
        trail[g_data->trail_count++] = next_var;
        continue;

handle_conflict:
        g_data->conflicts_found++;
        if (g_data->decision_level_current == 0) return; // UNSATISFIABLE

        // 3. BACKTRACK: Undo assignments until we find the decision to flip.
        while (g_data->trail_count > 0) {
            int var_to_undo = trail[g_data->trail_count - 1];
            g_data->trail_count--;
            
            // Check if this was the decision variable at the current conflict level.
            // Our simple strategy always decides TRUE (1), so we flip it to FALSE (-1).
            if (decision_level[var_to_undo] == g_data->decision_level_current) {
                assignment[var_to_undo] = -1; // Flip the decision
                g_data->decision_level_current--;
                decision_level[var_to_undo] = g_data->decision_level_current;
                trail[g_data->trail_count++] = var_to_undo;
                propagation_head = g_data->trail_count - 1;
                goto next_iteration;
            }
            
            assignment[var_to_undo] = 0;
            decision_level[var_to_undo] = 0;
        }
    next_iteration:;
    }
}

void cleanup() {
    if (g_data) {
        if (g_data->clauses) {
            for (int i = 0; i < g_data->num_clauses; ++i) {
                free(g_data->clauses[i].literals);
            }
            free(g_data->clauses);
        }
        free(g_data->assignment);
        free(g_data->decision_level);
        free(g_data->trail);
        free(g_data);
    }
}
