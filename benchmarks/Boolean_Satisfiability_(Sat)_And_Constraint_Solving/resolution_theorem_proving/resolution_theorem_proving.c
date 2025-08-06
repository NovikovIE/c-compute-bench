/**
 * @file resolution_theorem_proving.c
 * @brief A benchmark simulating a simplified resolution-based theorem prover for SAT.
 *
 * Description:
 * This program implements a core component of many constraint solvers: the
 * resolution rule. Given a set of boolean clauses (a formula in Conjunctive
 * Normal Form), the program repeatedly applies the resolution rule to generate
 * new clauses. The goal in a real solver is to derive an empty clause, which
 * proves the original formula is unsatisfiable (UNSAT).
 *
 * The resolution rule states that from two clauses (A v L) and (B v ~L),
 * we can infer a new clause (A v B), where A and B are disjunctions of other
 * literals and L is a literal. This process continues until no new clauses can
 * be derived (saturation) or an empty clause is found.
 *
 * The benchmark's computational workload comes from:
 * 1. The O(N^2) process of pairing up existing clauses to find resolution opportunities.
 * 2. The creation, sorting, and simplification of new clauses (resolvents).
 * 3. The O(N) check to ensure a newly derived clause is unique before adding it.
 * 4. Dynamic memory management as the set of clauses grows.
 *
 * This workload is representative of many NP-complete problems in logic and AI.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

// --- MERSENNE TWISTER (DO NOT MODIFY) ---
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

// --- BENCHMARK DATA STRUCTURES AND GLOBALS ---

// A literal is an integer. Positive for a variable, negative for its negation.
// e.g., variable 5 is literal 5, NOT(variable 5) is literal -5.
typedef struct {
    int* literals;
    int num_literals;
} Clause;

// Global parameters from command line
int NUM_INITIAL_CLAUSES;
int NUM_VARIABLES;
int MAX_LITERALS_PER_CLAUSE;
int MAX_TOTAL_CLAUSES;

// Global data structures for the knowledge base
Clause **clauses;
int num_clauses;
int clause_capacity;

// Benchmark result accumulator
int final_result;

// --- HELPER FUNCTIONS ---

// qsort comparison for integers
int compare_ints(const void* a, const void* b) {
    int arg1 = *(const int*)a;
    int arg2 = *(const int*)b;
    if (arg1 < arg2) return -1;
    if (arg1 > arg2) return 1;
    return 0;
}

// Allocates memory for a clause and its literals array
Clause* create_clause(int num_lits) {
    Clause* c = (Clause*)malloc(sizeof(Clause));
    if (!c) { perror("Failed to allocate clause"); exit(1); }
    c->num_literals = num_lits;
    if (num_lits > 0) {
        c->literals = (int*)malloc(sizeof(int) * num_lits);
        if (!c->literals) { perror("Failed to allocate literals"); exit(1); }
    } else {
        c->literals = NULL;
    }
    return c;
}

// Frees a clause and its internal literals array
void free_clause(Clause* c) {
    if (c) {
        free(c->literals);
        free(c);
    }
}

// Adds a new clause to the global knowledge base, resizing if necessary
void add_clause_to_kb(Clause* c) {
    if (num_clauses >= clause_capacity) {
        clause_capacity *= 2;
        clauses = (Clause**)realloc(clauses, sizeof(Clause*) * clause_capacity);
        if (!clauses) { perror("Failed to reallocate clauses array"); exit(1); }
    }
    clauses[num_clauses++] = c;
}

// Checks if two clauses are logically equivalent.
// Assumes literals within each clause are sorted.
int clauses_are_equal(Clause* c1, Clause* c2) {
    if (c1->num_literals != c2->num_literals) {
        return 0;
    }
    if (c1->num_literals == 0) return 1; // Both are empty clauses
    // Clauses are sorted, so a direct memory comparison is sufficient.
    return memcmp(c1->literals, c2->literals, c1->num_literals * sizeof(int)) == 0;
}

// Checks if a given clause already exists in the global knowledge base.
int is_clause_in_kb(Clause* new_clause) {
    for (int i = 0; i < num_clauses; ++i) {
        if (clauses_are_equal(clauses[i], new_clause)) {
            return 1;
        }
    }
    return 0;
}

// Tries to resolve two clauses, c1 and c2.
// Returns a new clause (the resolvent) if successful, otherwise NULL.
Clause* resolve(Clause* c1, Clause* c2) {
    for (int i = 0; i < c1->num_literals; ++i) {
        for (int j = 0; j < c2->num_literals; ++j) {
            // Check for a pivot: a literal in c1 and its negation in c2
            if (c1->literals[i] == -c2->literals[j]) {
                int potential_size = (c1->num_literals - 1) + (c2->num_literals - 1);
                if (potential_size == 0) {
                    return create_clause(0); // Empty clause found
                }

                int* temp_buffer = (int*)malloc(sizeof(int) * potential_size);
                if (!temp_buffer) { perror("Failed to allocate temp buffer"); exit(1); }
                int current_size = 0;

                // Copy literals from c1, skipping the pivot
                for (int k = 0; k < c1->num_literals; ++k) {
                    if (k != i) temp_buffer[current_size++] = c1->literals[k];
                }
                // Copy literals from c2, skipping the pivot
                for (int k = 0; k < c2->num_literals; ++k) {
                    if (k != j) temp_buffer[current_size++] = c2->literals[k];
                }
                
                qsort(temp_buffer, current_size, sizeof(int), compare_ints);

                // Check for tautology (e.g., A v ~A) and remove duplicates in a single pass
                if (current_size > 0) {
                    int unique_idx = 0;
                    for (int k = 1; k < current_size; ++k) {
                        if (temp_buffer[k] == -temp_buffer[unique_idx]) { // Tautology
                            free(temp_buffer);
                            return NULL;
                        }
                        if (temp_buffer[k] != temp_buffer[unique_idx]) {
                            temp_buffer[++unique_idx] = temp_buffer[k];
                        }
                    }
                    current_size = unique_idx + 1;
                }

                Clause* resolvent = create_clause(current_size);
                memcpy(resolvent->literals, temp_buffer, current_size * sizeof(int));
                free(temp_buffer);
                return resolvent;
            }
        }
    }
    return NULL; // No resolving literal pair found
}


// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_initial_clauses> <num_variables> <max_literals_per_clause> <seed>\n", argv[0]);
        exit(1);
    }

    NUM_INITIAL_CLAUSES = atoi(argv[1]);
    NUM_VARIABLES = atoi(argv[2]);
    MAX_LITERALS_PER_CLAUSE = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    
    // Safety cap to prevent runaway execution time and memory usage
    MAX_TOTAL_CLAUSES = NUM_INITIAL_CLAUSES * 10; 

    mt_seed(seed);

    clause_capacity = NUM_INITIAL_CLAUSES > 0 ? NUM_INITIAL_CLAUSES * 2 : 10;
    clauses = (Clause**)malloc(sizeof(Clause*) * clause_capacity);
    if (!clauses) { perror("Failed to allocate initial clauses array"); exit(1); }
    num_clauses = 0;
    
    char* used_vars = (char*)malloc(sizeof(char) * (NUM_VARIABLES + 1));
    if (!used_vars) { perror("Failed to allocate used_vars"); exit(1); };

    for (int i = 0; i < NUM_INITIAL_CLAUSES; ++i) {
        int num_lits = 1 + (mt_rand() % MAX_LITERALS_PER_CLAUSE);
        if (num_lits > NUM_VARIABLES) num_lits = NUM_VARIABLES;
        
        Clause* new_clause = create_clause(num_lits);
        memset(used_vars, 0, sizeof(char) * (NUM_VARIABLES + 1));
        
        for (int j = 0; j < num_lits; ++j) {
            int var;
            do {
                var = 1 + (mt_rand() % NUM_VARIABLES);
            } while (used_vars[var]);
            
            used_vars[var] = 1;
            int sign = (mt_rand() % 2 == 0) ? 1 : -1;
            new_clause->literals[j] = var * sign;
        }

        qsort(new_clause->literals, new_clause->num_literals, sizeof(int), compare_ints);
        add_clause_to_kb(new_clause);
    }
    
    free(used_vars);
    final_result = 0;
}

void run_computation() {
    int newly_added_in_round;
    int keep_running = 1;

    do {
        newly_added_in_round = 0;
        int clauses_at_start_of_round = num_clauses;

        for (int i = 0; i < clauses_at_start_of_round && keep_running; ++i) {
            for (int j = i + 1; j < clauses_at_start_of_round; ++j) {
                
                Clause* resolvent = resolve(clauses[i], clauses[j]); 
                if (resolvent) {
                    if (!is_clause_in_kb(resolvent)) {
                        add_clause_to_kb(resolvent);
                        newly_added_in_round++;
                        if (resolvent->num_literals == 0 || num_clauses >= MAX_TOTAL_CLAUSES) {
                            keep_running = 0; 
                            break; // Exit inner loop
                        }
                    } else {
                        free_clause(resolvent);
                    }
                }
            }
        }

    } while (newly_added_in_round > 0 && keep_running);

    final_result = num_clauses - NUM_INITIAL_CLAUSES;
}

void cleanup() {
    for (int i = 0; i < num_clauses; ++i) {
        free_clause(clauses[i]);
    }
    free(clauses);
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

    // Print result to stdout
    printf("%d\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
