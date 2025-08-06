#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- DO NOT MODIFY: Mersenne Twister Generator ---
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
// --- END Mersenne Twister ---

// --- Data Structures ---
enum TermType { CONSTANT, VARIABLE, FUNCTION };

typedef struct Term {
    enum TermType type;
    int id;
    int arity;
    struct Term** args;
} Term;

typedef struct {
    int predicate_id;
    _Bool is_negated;
    int arity;
    Term** args;
} Literal;

typedef struct {
    Literal** literals;
    int num_literals;
} Clause;

// --- Global Benchmark Parameters & Data ---
int num_initial_clauses;
int term_depth;
int literals_per_clause;
int resolution_steps;

Clause** clause_set = NULL;
int clause_count = 0;

volatile int final_result = 0;

// --- Function Prototypes ---
void setup_benchmark(int argc, char *argv[]);
void run_computation();
void cleanup();
Term* create_random_term(int depth);
void free_term(Term* t);
int simulate_unify_terms(Term* t1, Term* t2);

// --- Benchmark Implementation ---

Term* create_random_term(int depth) {
    Term* term = (Term*)malloc(sizeof(Term));
    if (!term) {
        perror("Failed to allocate term");
        exit(EXIT_FAILURE);
    }

    if (depth <= 0 || (mt_rand() % 4 == 0 && depth < term_depth)) {
        term->type = (mt_rand() % 2 == 0) ? CONSTANT : VARIABLE;
        term->id = mt_rand() % 1000;
        term->arity = 0;
        term->args = NULL;
    } else {
        term->type = FUNCTION;
        term->id = mt_rand() % 50;
        term->arity = 1 + (mt_rand() % 3);
        term->args = (Term**)malloc(term->arity * sizeof(Term*));
        if (!term->args) {
            perror("Failed to allocate term args");
            exit(EXIT_FAILURE);
        }
        for (int i = 0; i < term->arity; i++) {
            term->args[i] = create_random_term(depth - 1);
        }
    }
    return term;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <num_initial_clauses> <term_depth> <literals_per_clause> <resolution_steps> <seed>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    num_initial_clauses = atoi(argv[1]);
    term_depth = atoi(argv[2]);
    literals_per_clause = atoi(argv[3]);
    resolution_steps = atoi(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    mt_seed(seed);

    clause_set = (Clause**)malloc(num_initial_clauses * sizeof(Clause*));
    if (!clause_set) {
        perror("Failed to allocate clause set");
        exit(EXIT_FAILURE);
    }

    clause_count = num_initial_clauses;
    for (int i = 0; i < clause_count; i++) {
        clause_set[i] = (Clause*)malloc(sizeof(Clause));
        if (!clause_set[i]) {
            perror("Failed to allocate clause");
            exit(EXIT_FAILURE);
        }
        
        clause_set[i]->num_literals = literals_per_clause;
        clause_set[i]->literals = (Literal**)malloc(literals_per_clause * sizeof(Literal*));
        if (!clause_set[i]->literals) {
            perror("Failed to allocate literals array");
            exit(EXIT_FAILURE);
        }

        for (int j = 0; j < literals_per_clause; j++) {
            Literal* lit = (Literal*)malloc(sizeof(Literal));
            if (!lit) {
                perror("Failed to allocate literal");
                exit(EXIT_FAILURE);
            }
            lit->predicate_id = mt_rand() % 100;
            lit->is_negated = mt_rand() % 2;
            lit->arity = 1 + (mt_rand() % 3);
            lit->args = (Term**)malloc(lit->arity * sizeof(Term*));
             if (!lit->args) {
                perror("Failed to allocate literal args");
                exit(EXIT_FAILURE);
            }
            for (int k = 0; k < lit->arity; k++) {
                lit->args[k] = create_random_term(term_depth);
            }
            clause_set[i]->literals[j] = lit;
        }
    }
    final_result = 0;
}

int simulate_unify_terms(Term* t1, Term* t2) {
    if (!t1 || !t2) return 0;
    
    int cost = 1;

    if (t1->type == VARIABLE || t2->type == VARIABLE) {
        return cost;
    }

    if (t1->type != t2->type || t1->id != t2->id || t1->arity != t2->arity) {
        return 0; // Mismatch, unification fails
    }

    for (int i = 0; i < t1->arity; i++) {
        int sub_cost = simulate_unify_terms(t1->args[i], t2->args[i]);
        if (sub_cost == 0) {
            return 0; // A nested term failed to unify
        }
        cost += sub_cost;
    }
    
    return cost;
}

void run_computation() {
    int result = 0;
    for (int step = 0; step < resolution_steps; step++) {
        int idx1 = mt_rand() % clause_count;
        int idx2 = mt_rand() % clause_count;
        if (idx1 == idx2) continue;

        Clause* c1 = clause_set[idx1];
        Clause* c2 = clause_set[idx2];
        _Bool found_resolution = 0;

        for (int i = 0; i < c1->num_literals && !found_resolution; i++) {
            for (int j = 0; j < c2->num_literals && !found_resolution; j++) {
                Literal* l1 = c1->literals[i];
                Literal* l2 = c2->literals[j];

                if (l1->predicate_id == l2->predicate_id && l1->is_negated != l2->is_negated) {
                    if (l1->arity == l2->arity) {
                        int unification_cost = 0;
                        _Bool can_unify = 1;
                        for (int k = 0; k < l1->arity; k++) {
                           int sub_cost = simulate_unify_terms(l1->args[k], l2->args[k]);
                           if (sub_cost == 0) {
                               can_unify = 0;
                               break;
                           }
                           unification_cost += sub_cost;
                        }

                        if (can_unify) {
                            result += unification_cost;
                            result += (c1->num_literals - 1) + (c2->num_literals - 1);
                            found_resolution = 1;
                        }
                    }
                }
            }
        }
    }
    final_result = result;
}


void free_term(Term* t) {
    if (!t) return;
    if (t->type == FUNCTION) {
        for (int i = 0; i < t->arity; i++) {
            free_term(t->args[i]);
        }
        if (t->args) free(t->args);
    }
    free(t);
}

void cleanup() {
    if (!clause_set) return;

    for (int i = 0; i < clause_count; i++) {
        Clause* c = clause_set[i];
        if (!c) continue;
        
        for (int j = 0; j < c->num_literals; j++) {
            Literal* lit = c->literals[j];
            if (!lit) continue;
            
            for (int k = 0; k < lit->arity; k++) {
                free_term(lit->args[k]);
            }
            if (lit->args) free(lit->args);
            free(lit);
        }
        if (c->literals) free(c->literals);
        free(c);
    }
    free(clause_set);
}

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
