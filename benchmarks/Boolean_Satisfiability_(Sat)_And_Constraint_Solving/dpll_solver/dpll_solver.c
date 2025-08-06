#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>


// --- Mersenne Twister (Do Not Modify - Include This Verbatim) ---
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

// --- Benchmark Globals & Data Structures ---
typedef struct {
    int num_variables;
    int num_clauses;
    int literals_per_clause;
    int** clauses; // Represents the CNF formula
    int* assignment;  // Stores variable assignments
    int result;       // 1 for SAT, 0 for UNSAT
} BenchmarkData;

static BenchmarkData g_data;

// Enum for variable states
enum { VAR_UNASSIGNED = -1, VAR_FALSE = 0, VAR_TRUE = 1 };

// Forward declaration for the recursive solver
bool dpll_solve();

// --- Data Setup --- 
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_variables> <num_clauses> <clause_to_variable_ratio> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_variables = atoi(argv[1]);
    g_data.num_clauses = atoi(argv[2]);
    g_data.literals_per_clause = atoi(argv[3]);
    uint32_t seed = (uint32_t)strtoul(argv[4], NULL, 10);
    
    mt_seed(seed);

    if (g_data.literals_per_clause > g_data.num_variables) {
        fprintf(stderr, "FATAL: literals_per_clause cannot be greater than num_variables.\n");
        exit(1);
    }

    // Allocate memory for clauses (an array of int arrays)
    g_data.clauses = (int**)malloc(g_data.num_clauses * sizeof(int*));
    if (g_data.clauses == NULL) { exit(1); }

    // Generate clauses
    int* var_indices = (int*)malloc(g_data.num_variables * sizeof(int));
    if (var_indices == NULL) { exit(1); }
    for (int i = 0; i < g_data.num_variables; i++) {
        var_indices[i] = i + 1;
    }

    for (int i = 0; i < g_data.num_clauses; i++) {
        g_data.clauses[i] = (int*)malloc(g_data.literals_per_clause * sizeof(int));
        if (g_data.clauses[i] == NULL) { exit(1); }

        // Fisher-Yates shuffle to pick unique variables for the clause
        for (int j = g_data.num_variables - 1; j > 0; j--) {
            int k = mt_rand() % (j + 1);
            int temp = var_indices[j];
            var_indices[j] = var_indices[k];
            var_indices[k] = temp;
        }

        for(int j = 0; j < g_data.literals_per_clause; j++) {
            int literal = var_indices[j];
            if (mt_rand() % 2) { // Randomly negate the literal
                literal = -literal;
            }
            g_data.clauses[i][j] = literal;
        }
    }
    free(var_indices);

    // Allocate and initialize assignment array
    // Index 0 is unused, variables are 1-based.
    g_data.assignment = (int*)malloc((g_data.num_variables + 1) * sizeof(int));
    if (g_data.assignment == NULL) { exit(1); }
    for (int i = 0; i <= g_data.num_variables; i++) {
        g_data.assignment[i] = VAR_UNASSIGNED;
    }
}

// --- Computation ---
void run_computation() {
    g_data.result = dpll_solve() ? 1 : 0;
}

// --- Memory Cleanup ---
void cleanup() {
    for (int i = 0; i < g_data.num_clauses; i++) {
        free(g_data.clauses[i]);
    }
    free(g_data.clauses);
    free(g_data.assignment);
}

// --- Core DPLL Algorithm ---
bool dpll_solve() {
    // --- 1. Unit Propagation ---
    // This simplified version of DPLL integrates unit propagation into the main check.
    // A more optimized version would have a separate, iterative unit propagation loop.

    // --- 2. Check formula status (SAT, UNSAT, or Undetermined) ---
    bool all_clauses_satisfied = true;
    for (int i = 0; i < g_data.num_clauses; i++) {
        bool clause_satisfied = false;
        bool all_literals_false = true;
        for (int j = 0; j < g_data.literals_per_clause; j++) {
            int literal = g_data.clauses[i][j];
            int var = abs(literal);
            int value = g_data.assignment[var];

            if (value == VAR_UNASSIGNED) {
                all_literals_false = false;
            } else if ((literal > 0 && value == VAR_TRUE) || (literal < 0 && value == VAR_FALSE)) {
                clause_satisfied = true;
                all_literals_false = false;
                break; // Clause is satisfied
            }
        }

        if (!clause_satisfied) {
            all_clauses_satisfied = false;
            if (all_literals_false) {
                return false; // Conflict found
            }
        }
    }

    if (all_clauses_satisfied) {
        return true; // Solution found
    }

    // --- 3. Branching ---
    // Find the first unassigned variable to branch on
    int branch_var = 0;
    for (int i = 1; i <= g_data.num_variables; i++) {
        if (g_data.assignment[i] == VAR_UNASSIGNED) {
            branch_var = i;
            break;
        }
    }

    if (branch_var == 0) {
        // This case should be covered by checks above but is a safeguard
        return false;
    }

    // Store assignment state for backtracking
    size_t assignment_size = (g_data.num_variables + 1) * sizeof(int);
    int* assignment_backup = (int*)malloc(assignment_size);
    if (assignment_backup == NULL) { exit(1); }
    memcpy(assignment_backup, g_data.assignment, assignment_size);

    // Branch 1: Try assigning variable to TRUE
    g_data.assignment[branch_var] = VAR_TRUE;
    if (dpll_solve()) {
        free(assignment_backup);
        return true;
    }
    
    // Backtrack from TRUE branch: restore assignment
    memcpy(g_data.assignment, assignment_backup, assignment_size);

    // Branch 2: Try assigning variable to FALSE
    g_data.assignment[branch_var] = VAR_FALSE;
    if (dpll_solve()) {
        free(assignment_backup);
        return true;
    }

    // Backtrack from FALSE branch: restore assignment
    memcpy(g_data.assignment, assignment_backup, assignment_size);
    free(assignment_backup);

    return false; // Both branches failed
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

    // Print result to stdout
    printf("%d\n", g_data.result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
