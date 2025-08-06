#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <string.h>

// --- Mersenne Twister (MT19937) Generator ---
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

// --- Global Benchmark Data ---
int num_variables;
int num_clauses;
int num_hard_clauses;

// Each clause is an array of literals (integers), terminated by 0.
// A positive integer 'v' means variable v-1 is true.
// A negative integer '-v' means variable v-1 is false.
int** clauses;

// Current assignment for variables: 0 for false, 1 for true.
int* assignment;

// The computational result.
int final_result;

// --- Forward Declarations ---
void setup_benchmark(int argc, char *argv[]);
void run_computation();
void cleanup();
void backtracking_search(int k);

// --- Main --- 
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


// --- Function Implementations ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_variables> <num_clauses> <num_hard_clauses> <seed>\n", argv[0]);
        exit(1);
    }

    num_variables = atoi(argv[1]);
    num_clauses = atoi(argv[2]);
    num_hard_clauses = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    
    if (num_hard_clauses > num_clauses) {
        fprintf(stderr, "Error: num_hard_clauses cannot be greater than num_clauses.\n");
        exit(1);
    }

    mt_seed(seed);

    clauses = (int**)malloc(sizeof(int*) * num_clauses);
    assignment = (int*)malloc(sizeof(int) * num_variables);
    int* vars_in_clause = (int*)malloc(sizeof(int) * num_variables);

    for (int i = 0; i < num_clauses; ++i) {
        // Clauses have 2 to 4 literals to keep them non-trivial.
        int clause_len = 2 + (mt_rand() % 3);
        clauses[i] = (int*)malloc(sizeof(int) * (clause_len + 1));
        memset(vars_in_clause, 0, sizeof(int) * num_variables);
        
        int current_len = 0;
        while(current_len < clause_len) {
            int var_idx = mt_rand() % num_variables;
            if (vars_in_clause[var_idx] == 0) {
                int literal = var_idx + 1;
                if (mt_rand() % 2) { // Randomly negate
                    literal = -literal;
                }
                clauses[i][current_len++] = literal;
                vars_in_clause[var_idx] = 1;
            }
        }
        clauses[i][current_len] = 0; // Null-terminate the clause
    }

    free(vars_in_clause);
}

void backtracking_search(int k) {
    // Base Case: all variables have been assigned.
    if (k == num_variables) {
        // Step 1: Check if the current assignment satisfies all hard clauses.
        for (int i = 0; i < num_hard_clauses; ++i) {
            int is_satisfied = 0;
            for (int* lit_ptr = clauses[i]; *lit_ptr != 0; ++lit_ptr) {
                int literal = *lit_ptr;
                int var_idx = abs(literal) - 1;
                int is_positive = literal > 0;
                if ((is_positive && assignment[var_idx] == 1) || 
                    (!is_positive && assignment[var_idx] == 0)) {
                    is_satisfied = 1;
                    break;
                }
            }
            if (!is_satisfied) {
                return; // Hard clause violated, this assignment is not valid.
            }
        }

        // Step 2: If all hard clauses are satisfied, count satisfied soft clauses.
        int current_soft_satisfied = 0;
        for (int i = num_hard_clauses; i < num_clauses; ++i) {
             int is_satisfied = 0;
            for (int* lit_ptr = clauses[i]; *lit_ptr != 0; ++lit_ptr) {
                int literal = *lit_ptr;
                int var_idx = abs(literal) - 1;
                int is_positive = literal > 0;
                if ((is_positive && assignment[var_idx] == 1) || 
                    (!is_positive && assignment[var_idx] == 0)) {
                    is_satisfied = 1;
                    break;
                }
            }
            if (is_satisfied) {
                current_soft_satisfied++;
            }
        }

        // Step 3: Update the global best result.
        if (final_result < current_soft_satisfied) {
            final_result = current_soft_satisfied;
        }
        return;
    }

    // Recursive Step: Branch on the current variable 'k'.
    // Branch 1: Try assigning variable 'k' to true (1).
    assignment[k] = 1;
    backtracking_search(k + 1);

    // Branch 2: Try assigning variable 'k' to false (0).
    assignment[k] = 0;
    backtracking_search(k + 1);
}

void run_computation() {
    // Initialize result to -1, indicating no valid solution has been found yet.
    final_result = -1; 
    backtracking_search(0);
}

void cleanup() {
    for (int i = 0; i < num_clauses; ++i) {
        free(clauses[i]);
    }
    free(clauses);
    free(assignment);
}
