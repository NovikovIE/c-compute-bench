#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <stdbool.h>

// --- MERSENNE TWISTER (Do Not Modify - Include This Verbatim) ---
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

// --- BENCHMARK DATA STRUCTURES ---

typedef enum { VAR, FUNC } TermType;

// Represents a logical term in first-order logic.
typedef struct Term {
    TermType type;
    union {
        // A variable, identified by an integer ID.
        int var_id;
        // A function symbol with arguments.
        struct {
            int func_id;    // Identifier for the function symbol (e.g., f, g, a, b).
            int arity;      // Number of arguments.
            struct Term** args; // Array of pointers to argument terms.
        } func;
    } data;
} Term;

// Global structure to hold all benchmark data and parameters.
typedef struct {
    int num_terms;
    int max_term_depth;
    int num_variables;

    // Pairs of terms to be unified.
    Term **term_pairs1;
    Term **term_pairs2;

    // Final result accumulator.
    int successful_unifications;
} BenchmarkData;

static BenchmarkData g_data;

// --- FORWARD DECLARATIONS ---
Term* generate_term_recursive(int current_depth);
void free_term(Term* t);
bool unify(Term* t1, Term* t2, Term** substitution);
bool occurs_check(int var_id, Term* t, Term** substitution);

// --- BENCHMARK IMPLEMENTATION ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_terms> <max_term_depth> <num_variables> <seed>\n", argv[0]);
        exit(1);
    }
    g_data.num_terms = atoi(argv[1]);
    g_data.max_term_depth = atoi(argv[2]);
    g_data.num_variables = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    mt_seed(seed);

    g_data.term_pairs1 = (Term**)malloc(g_data.num_terms * sizeof(Term*));
    g_data.term_pairs2 = (Term**)malloc(g_data.num_terms * sizeof(Term*));
    if (!g_data.term_pairs1 || !g_data.term_pairs2) {
        fprintf(stderr, "FATAL: Memory allocation failed for term pairs.\n");
        exit(1);
    }

    for (int i = 0; i < g_data.num_terms; i++) {
        g_data.term_pairs1[i] = generate_term_recursive(0);
        g_data.term_pairs2[i] = generate_term_recursive(0);
    }

    g_data.successful_unifications = 0;
}

void run_computation() {
    int success_count = 0;
    Term** substitution = (Term**)malloc(g_data.num_variables * sizeof(Term*));
    if (!substitution) {
        fprintf(stderr, "FATAL: Memory allocation failed for substitution array.\n");
        exit(1);
    }

    for (int i = 0; i < g_data.num_terms; i++) {
        // Reset substitution for each new unification problem.
        memset(substitution, 0, g_data.num_variables * sizeof(Term*));
        
        if (unify(g_data.term_pairs1[i], g_data.term_pairs2[i], substitution)) {
            success_count++;
        }
    }

    free(substitution);
    g_data.successful_unifications = success_count;
}

void cleanup() {
    for (int i = 0; i < g_data.num_terms; i++) {
        free_term(g_data.term_pairs1[i]);
        free_term(g_data.term_pairs2[i]);
    }
    free(g_data.term_pairs1);
    free(g_data.term_pairs2);
}

// --- HELPER FUNCTIONS ---

// Recursively generates a random term.
Term* generate_term_recursive(int current_depth) {
    Term* t = (Term*)malloc(sizeof(Term));
    if (!t) {
        fprintf(stderr, "FATAL: Memory allocation failed for a term node.\n");
        exit(1);
    }

    // Bias towards terminating with a variable as depth increases.
    bool terminate = (current_depth >= g_data.max_term_depth) || ((mt_rand() % 100) < (25 + current_depth * 5));

    if (terminate || g_data.num_variables == 0) {
        t->type = VAR;
        t->data.var_id = mt_rand() % g_data.num_variables;
    } else {
        t->type = FUNC;
        t->data.func.func_id = mt_rand() % 50; // 50 different function symbols
        t->data.func.arity = 1 + (mt_rand() % 3); // Arity from 1 to 3
        t->data.func.args = (Term**)malloc(t->data.func.arity * sizeof(Term*));
        if (!t->data.func.args) {
            fprintf(stderr, "FATAL: Memory allocation failed for term arguments.\n");
            exit(1);
        }
        for (int i = 0; i < t->data.func.arity; i++) {
            t->data.func.args[i] = generate_term_recursive(current_depth + 1);
        }
    }
    return t;
}

// Recursively frees memory allocated for a term.
void free_term(Term* t) {
    if (t == NULL) return;
    if (t->type == FUNC) {
        for (int i = 0; i < t->data.func.arity; i++) {
            free_term(t->data.func.args[i]);
        }
        free(t->data.func.args);
    }
    free(t);
}

// Prevents unifying a variable with a term containing that same variable (e.g., X = f(X)).
bool occurs_check(int var_id, Term* t, Term** substitution) {
    // Dereference t to its ultimate binding
    while(t->type == VAR && substitution[t->data.var_id] != NULL) {
        t = substitution[t->data.var_id];
    }
    
    if (t->type == VAR) {
        return t->data.var_id == var_id;
    }

    if (t->type == FUNC) {
        for (int i = 0; i < t->data.func.arity; i++) {
            if (occurs_check(var_id, t->data.func.args[i], substitution)) {
                return true;
            }
        }
    }
    return false;
}

// Core unification algorithm.
bool unify(Term* t1, Term* t2, Term** substitution) {
    // Dereference both terms to their final values.
    while (t1->type == VAR && substitution[t1->data.var_id] != NULL) {
        t1 = substitution[t1->data.var_id];
    }
    while (t2->type == VAR && substitution[t2->data.var_id] != NULL) {
        t2 = substitution[t2->data.var_id];
    }

    if (t1 == t2) return true;

    if (t1->type == VAR) {
        if (occurs_check(t1->data.var_id, t2, substitution)) return false;
        substitution[t1->data.var_id] = t2;
        return true;
    }

    if (t2->type == VAR) {
        if (occurs_check(t2->data.var_id, t1, substitution)) return false;
        substitution[t2->data.var_id] = t1;
        return true;
    }

    if (t1->type == FUNC && t2->type == FUNC) {
        if (t1->data.func.func_id != t2->data.func.func_id || t1->data.func.arity != t2->data.func.arity) {
            return false; // Mismatch in function symbol or arity.
        }
        for (int i = 0; i < t1->data.func.arity; i++) {
            if (!unify(t1->data.func.args[i], t2->data.func.args[i], substitution)) {
                return false; // Unification of arguments failed.
            }
        }
        return true;
    }
    
    return false; // Should not be reached.
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
    printf("%d\n", g_data.successful_unifications);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
