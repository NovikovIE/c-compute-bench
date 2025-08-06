#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <string.h>
#include <limits.h>
#include <math.h>

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

// --- BENCHMARK DATA AND FUNCTIONS ---

#define CLAUSE_SIZE 3 // Using 3-SAT problems

// Global struct to hold all benchmark data
typedef struct {
    int num_variables;
    int num_clauses;
    long max_flips;
    int max_tries;

    // CNF formula: clauses[i][j] is the j-th literal of the i-th clause.
    // A literal is a non-zero integer. Positive for a variable, negative for its negation.
    // Variable `v` (0-indexed) is represented by integer `v+1`.
    int** clauses;
    bool* assignment; // Current truth assignment for variables
    
    // Final result: number of times a satisfying assignment was found
    int found_solutions_count;
} BenchmarkState;

static BenchmarkState g_data;

// --- Function Prototypes ---
void setup_benchmark(int argc, char* argv[]);
void run_computation();
void cleanup();

// --- Data Setup ---
void setup_benchmark(int argc, char* argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s num_variables num_clauses max_flips max_tries seed\n", argv[0]);
        exit(1);
    }

    g_data.num_variables = atoi(argv[1]);
    g_data.num_clauses = atoi(argv[2]);
    g_data.max_flips = atol(argv[3]);
    g_data.max_tries = atoi(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    mt_seed(seed);

    // Allocate memory for the clauses and the assignment
    g_data.clauses = (int**)malloc(g_data.num_clauses * sizeof(int*));
    if (!g_data.clauses) { perror("malloc clauses"); exit(1); }
    for (int i = 0; i < g_data.num_clauses; ++i) {
        g_data.clauses[i] = (int*)malloc(CLAUSE_SIZE * sizeof(int));
        if (!g_data.clauses[i]) { perror("malloc clause row"); exit(1); }
    }
    g_data.assignment = (bool*)malloc(g_data.num_variables * sizeof(bool));
    if (!g_data.assignment) { perror("malloc assignment"); exit(1); }

    // Generate a random 3-SAT instance
    for (int i = 0; i < g_data.num_clauses; ++i) {
        int vars[CLAUSE_SIZE];
        for (int j = 0; j < CLAUSE_SIZE; ++j) {
            int var;
            bool unique;
            do {
                unique = true;
                var = (mt_rand() % g_data.num_variables) + 1;
                for (int k = 0; k < j; ++k) {
                    if (vars[k] == var) {
                        unique = false;
                        break;
                    }
                }
            } while (!unique);
            vars[j] = var;
            int sign = (mt_rand() % 2 == 0) ? 1 : -1;
            g_data.clauses[i][j] = var * sign;
        }
    }

    g_data.found_solutions_count = 0;
}

// --- Core Computation ---

// Checks if a clause is satisfied by the current global assignment.
static inline bool is_clause_satisfied(int clause_idx) {
    for (int i = 0; i < CLAUSE_SIZE; ++i) {
        int literal = g_data.clauses[clause_idx][i];
        int var_idx = abs(literal) - 1;
        bool is_negated = literal < 0;
        if (g_data.assignment[var_idx] != is_negated) {
            return true; // Literal is true, so clause is satisfied
        }
    }
    return false;
}

void run_computation() {
    int* unsatisfied_clauses = (int*)malloc(g_data.num_clauses * sizeof(int));
    if (!unsatisfied_clauses) { perror("malloc unsatisfied_clauses"); exit(1); }

    for (int try_num = 0; try_num < g_data.max_tries; ++try_num) {
        // Start with a random assignment
        for (int i = 0; i < g_data.num_variables; ++i) {
            g_data.assignment[i] = (mt_rand() % 2 == 0);
        }

        for (long flip_num = 0; flip_num < g_data.max_flips; ++flip_num) {
            // Find unsatisfied clauses
            int num_unsatisfied = 0;
            for (int i = 0; i < g_data.num_clauses; ++i) {
                if (!is_clause_satisfied(i)) {
                    unsatisfied_clauses[num_unsatisfied++] = i;
                }
            }

            if (num_unsatisfied == 0) {
                g_data.found_solutions_count++;
                break; // Found solution, start next try
            }

            // Pick a random unsatisfied clause
            int clause_to_fix_idx = unsatisfied_clauses[mt_rand() % num_unsatisfied];

            // Heuristic: with 50% probability, flip a random variable in the clause
            if (mt_rand() % 100 < 50) {
                int literal_to_flip = g_data.clauses[clause_to_fix_idx][mt_rand() % CLAUSE_SIZE];
                int var_idx_to_flip = abs(literal_to_flip) - 1;
                g_data.assignment[var_idx_to_flip] = !g_data.assignment[var_idx_to_flip];
            } else {
                // Greedy move: find the variable that breaks the fewest clauses when flipped
                int best_var_idx = -1;
                int min_break_count = INT_MAX;

                for (int i = 0; i < CLAUSE_SIZE; ++i) {
                    int var_idx = abs(g_data.clauses[clause_to_fix_idx][i]) - 1;
                    int current_break_count = 0;

                    // Flip temporarily to calculate break count
                    g_data.assignment[var_idx] = !g_data.assignment[var_idx];
                    
                    for (int c = 0; c < g_data.num_clauses; ++c) {
                        if (!is_clause_satisfied(c)) {
                           current_break_count++;
                        }
                    }

                    // Flip back
                    g_data.assignment[var_idx] = !g_data.assignment[var_idx];

                    if (current_break_count < min_break_count) {
                        min_break_count = current_break_count;
                        best_var_idx = var_idx;
                    }
                }
                g_data.assignment[best_var_idx] = !g_data.assignment[best_var_idx];
            }
        }
    }
    free(unsatisfied_clauses);
}

// --- Cleanup ---
void cleanup() {
    for (int i = 0; i < g_data.num_clauses; ++i) {
        free(g_data.clauses[i]);
    }
    free(g_data.clauses);
    free(g_data.assignment);
}

// --- Main Driver ---
int main(int argc, char* argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout
    printf("%d\n", g_data.found_solutions_count);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
