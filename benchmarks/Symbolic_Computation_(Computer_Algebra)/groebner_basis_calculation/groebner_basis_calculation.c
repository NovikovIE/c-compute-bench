/**
 * @file groebner_basis_calculation.c
 * @brief A benchmark simulating Groebner basis calculation using a simplified Buchberger's algorithm.
 *
 * This program represents a symbolic computation workload. The core of the benchmark
 * is the manipulation of multivariate polynomials. The complexity stems from the
 * generation of S-polynomials and their reduction, which can lead to "expression swell"—
 * a rapid increase in the number and complexity of intermediate polynomials.
 *
 * The benchmark creates an initial set of random polynomials and then iteratively
 * refines this set to form a Groebner basis. The primary operations are polynomial
 * arithmetic (multiplication, subtraction) over integer coefficients.
 *
 * Data representation:
 * - A Term is a coefficient and a list of integer exponents (e.g., 5 * x^2 * y^3 * z^0).
 * - A Polynomial is a sorted list of Terms.
 * - Lexicographical order is used for comparing monomials.
 *
 * The final result is the sum of all coefficients in the final computed basis,
 * which serves as a checksum to prevent dead-code elimination.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

// --- Mersenne Twister (MT19937) implementation ---
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
// --- End of MT19937 ---

// --- Benchmark Data Structures and Globals ---

// A term is coeff * x0^e0 * x1^e1 * ...
typedef struct {
    long long coefficient;
    int* exponents; // Array of size NUM_VARIABLES
} Term;

// A polynomial is a list of terms, sorted by monomial order
typedef struct {
    Term* terms;
    int num_terms;
    int capacity;
} Polynomial;

typedef struct {
    int p1_idx;
    int p2_idx;
} CriticalPair;

// Parameters
static int NUM_POLYNOMIALS;
static int NUM_VARIABLES;
static int MAX_DEGREE;

// Global state for the Groebner basis calculation
static Polynomial** G = NULL; // The basis, which will grow
static int g_size = 0;
static int g_capacity = 0;

static CriticalPair* critical_pairs = NULL;
static int pair_count = 0;
static int pair_capacity = 0;

// Final result to be printed to stdout
static long long final_result = 0;

// --- Forward declarations for polynomial operations ---
void free_poly(Polynomial* p);
void sort_poly_terms(Polynomial* p);
Polynomial* subtract_poly(const Polynomial* p1, const Polynomial* p2);
Polynomial* term_multiply_poly(const Polynomial* p, const Term* t);

// --- Utility functions ---

int compare_monomials(const int* exp1, const int* exp2) {
    for (int i = 0; i < NUM_VARIABLES; ++i) {
        if (exp1[i] > exp2[i]) return 1;
        if (exp1[i] < exp2[i]) return -1;
    }
    return 0;
}

int term_qsort_comparator(const void* a, const void* b) {
    const Term* t1 = (const Term*)a;
    const Term* t2 = (const Term*)b;
    return -compare_monomials(t1->exponents, t2->exponents); // Descending order
}

void add_term_to_poly(Polynomial* p, Term t) {
    if (p->num_terms >= p->capacity) {
        p->capacity = p->capacity == 0 ? 8 : p->capacity * 2;
        p->terms = (Term*)realloc(p->terms, p->capacity * sizeof(Term));
    }
    p->terms[p->num_terms++] = t;
}

void add_poly_to_basis(Polynomial* p) {
    if (g_size >= g_capacity) {
        g_capacity = g_capacity == 0 ? 8 : g_capacity * 2;
        G = (Polynomial**)realloc(G, g_capacity * sizeof(Polynomial*));
    }
    G[g_size++] = p;
}

void add_critical_pair(int i, int j) {
    if (pair_count >= pair_capacity) {
        pair_capacity = pair_capacity == 0 ? 16 : pair_capacity * 2;
        critical_pairs = (CriticalPair*)realloc(critical_pairs, pair_capacity * sizeof(CriticalPair));
    }
    critical_pairs[pair_count].p1_idx = i;
    critical_pairs[pair_count].p2_idx = j;
    pair_count++;
}

// --- Polynomial Arithmetic ---

void sort_poly_terms(Polynomial* p) {
    if (p->num_terms > 1) {
        qsort(p->terms, p->num_terms, sizeof(Term), term_qsort_comparator);
    }
}

void free_poly(Polynomial* p) {
    if (!p) return;
    for (int i = 0; i < p->num_terms; ++i) {
        free(p->terms[i].exponents);
    }
    free(p->terms);
    free(p);
}

Polynomial* create_poly() {
    Polynomial* p = (Polynomial*)calloc(1, sizeof(Polynomial));
    return p;
}

Polynomial* copy_poly(const Polynomial* p) {
    Polynomial* new_p = create_poly();
    for (int i = 0; i < p->num_terms; ++i) {
        Term new_t;
        new_t.coefficient = p->terms[i].coefficient;
        new_t.exponents = (int*)malloc(NUM_VARIABLES * sizeof(int));
        memcpy(new_t.exponents, p->terms[i].exponents, NUM_VARIABLES * sizeof(int));
        add_term_to_poly(new_p, new_t);
    }
    return new_p;
}

// Subtracts p2 from p1, returns a new polynomial.
// Assumes p1 and p2 have sorted terms.
Polynomial* subtract_poly(const Polynomial* p1, const Polynomial* p2) {
    Polynomial* res = create_poly();
    int i = 0, j = 0;

    while (i < p1->num_terms && j < p2->num_terms) {
        int cmp = compare_monomials(p1->terms[i].exponents, p2->terms[j].exponents);
        Term new_t;
        new_t.exponents = (int*)malloc(NUM_VARIABLES * sizeof(int));

        if (cmp > 0) {
            new_t.coefficient = p1->terms[i].coefficient;
            memcpy(new_t.exponents, p1->terms[i].exponents, NUM_VARIABLES * sizeof(int));
            i++;
        } else if (cmp < 0) {
            new_t.coefficient = -p2->terms[j].coefficient;
            memcpy(new_t.exponents, p2->terms[j].exponents, NUM_VARIABLES * sizeof(int));
            j++;
        } else {
            long long coeff = p1->terms[i].coefficient - p2->terms[j].coefficient;
            if (coeff != 0) {
                new_t.coefficient = coeff;
                memcpy(new_t.exponents, p1->terms[i].exponents, NUM_VARIABLES * sizeof(int));
            } else {
                free(new_t.exponents);
                i++; j++;
                continue;
            }
            i++; j++;
        }
        add_term_to_poly(res, new_t);
    }

    while (i < p1->num_terms) {
        Term new_t;
        new_t.exponents = (int*)malloc(NUM_VARIABLES * sizeof(int));
        new_t.coefficient = p1->terms[i].coefficient;
        memcpy(new_t.exponents, p1->terms[i].exponents, NUM_VARIABLES * sizeof(int));
        add_term_to_poly(res, new_t);
        i++;
    }
    while (j < p2->num_terms) {
        Term new_t;
        new_t.exponents = (int*)malloc(NUM_VARIABLES * sizeof(int));
        new_t.coefficient = -p2->terms[j].coefficient;
        memcpy(new_t.exponents, p2->terms[j].exponents, NUM_VARIABLES * sizeof(int));
        add_term_to_poly(res, new_t);
        j++;
    }
    return res;
}

Polynomial* term_multiply_poly(const Polynomial* p, const Term* t) {
    Polynomial* res = create_poly();
    for (int i = 0; i < p->num_terms; ++i) {
        Term new_t;
        new_t.coefficient = p->terms[i].coefficient * t->coefficient;
        new_t.exponents = (int*)malloc(NUM_VARIABLES * sizeof(int));
        for (int k = 0; k < NUM_VARIABLES; ++k) {
            new_t.exponents[k] = p->terms[i].exponents[k] + t->exponents[k];
        }
        add_term_to_poly(res, new_t);
    }
    return res;
}

// --- Core Algorithm ---

Polynomial* calculate_spoly(const Polynomial* p1, const Polynomial* p2) {
    if (p1->num_terms == 0 || p2->num_terms == 0) return create_poly();

    const Term* lt1 = &p1->terms[0];
    const Term* lt2 = &p2->terms[0];
    
    Term lcm_term1, lcm_term2;
    lcm_term1.exponents = (int*)malloc(NUM_VARIABLES * sizeof(int));
    lcm_term2.exponents = (int*)malloc(NUM_VARIABLES * sizeof(int));

    lcm_term1.coefficient = lt2->coefficient;
    lcm_term2.coefficient = lt1->coefficient;

    for (int i = 0; i < NUM_VARIABLES; ++i) {
        int exp_lcm = (lt1->exponents[i] > lt2->exponents[i]) ? lt1->exponents[i] : lt2->exponents[i];
        lcm_term1.exponents[i] = exp_lcm - lt1->exponents[i];
        lcm_term2.exponents[i] = exp_lcm - lt2->exponents[i];
    }

    Polynomial* s1 = term_multiply_poly(p1, &lcm_term1);
    Polynomial* s2 = term_multiply_poly(p2, &lcm_term2);
    Polynomial* spoly = subtract_poly(s1, s2);

    free(lcm_term1.exponents);
    free(lcm_term2.exponents);
    free_poly(s1);
    free_poly(s2);

    return spoly;
}

int is_monomial_divisible(const int* exp_num, const int* exp_den) {
    for (int i = 0; i < NUM_VARIABLES; ++i) {
        if (exp_num[i] < exp_den[i]) return 0;
    }
    return 1;
}

Polynomial* reduce_poly(Polynomial* p) {
    Polynomial* h = copy_poly(p);
    int reduction_occurred;
    do {
        reduction_occurred = 0;
        if (h->num_terms == 0) break;
        sort_poly_terms(h);
        
        const Term* lt_h = &h->terms[0];
        for (int i = 0; i < g_size; ++i) {
            if (G[i]->num_terms == 0) continue;
            const Term* lt_g = &G[i]->terms[0];

            if (is_monomial_divisible(lt_h->exponents, lt_g->exponents)) {
                Term quot_term;
                quot_term.exponents = (int*)malloc(NUM_VARIABLES * sizeof(int));
                quot_term.coefficient = lt_h->coefficient / lt_g->coefficient;
                // For integer arithmetic, we rely on having made polys monic, or this is simplified.
                // Here, we just assume it's divisible which might not be true for coefficients in general.
                // For this benchmark's integer simulation, we proceed as if it is.
                if (lt_g->coefficient != 0) {
                     quot_term.coefficient = lt_h->coefficient / lt_g->coefficient; 
                } else { // Should not happen with non-zero polys
                     quot_term.coefficient = 0;
                }

                for (int k = 0; k < NUM_VARIABLES; ++k) {
                    quot_term.exponents[k] = lt_h->exponents[k] - lt_g->exponents[k];
                }

                Polynomial* to_subtract = term_multiply_poly(G[i], &quot_term);
                Polynomial* next_h = subtract_poly(h, to_subtract);
                
                free_poly(h);
                free_poly(to_subtract);
                free(quot_term.exponents);
                h = next_h;
                reduction_occurred = 1;
                break; // Restart reduction with the new h
            }
        }
    } while (reduction_occurred);

    return h;
}


// --- Benchmark Setup, Run, Cleanup ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_polynomials num_variables max_degree seed\n", argv[0]);
        exit(1);
    }
    NUM_POLYNOMIALS = atoi(argv[1]);
    NUM_VARIABLES = atoi(argv[2]);
    MAX_DEGREE = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    mt_seed(seed);

    for (int i = 0; i < NUM_POLYNOMIALS; i++) {
        Polynomial* p = create_poly();
        int num_terms = (mt_rand() % (MAX_DEGREE * 2)) + 2;
        int total_degree_sum = 0;

        for (int j = 0; j < num_terms; j++) {
            Term t;
            t.coefficient = (mt_rand() % 10) - 5;
            if (t.coefficient == 0) t.coefficient = 1;

            t.exponents = (int*)calloc(NUM_VARIABLES, sizeof(int));
            int term_degree = 0;
            for (int k = 0; k < NUM_VARIABLES; k++) {
                if (term_degree < MAX_DEGREE) {
                    int exp = mt_rand() % (MAX_DEGREE - term_degree + 1);
                    t.exponents[k] = exp;
                    term_degree += exp;
                }
            }
            add_term_to_poly(p, t);
        }
        sort_poly_terms(p);
        // Make leading coefficient 1 or -1 to simplify reduction division
        if(p->num_terms > 0 && p->terms[0].coefficient != 1 && p->terms[0].coefficient != -1) {
            if (p->terms[0].coefficient > 1) p->terms[0].coefficient = 1;
            else if (p->terms[0].coefficient < -1) p->terms[0].coefficient = -1;
            // Not fully monic, but avoids division by large numbers
        }

        add_poly_to_basis(p);
    }

    for (int i = 0; i < g_size; i++) {
        for (int j = i + 1; j < g_size; j++) {
            add_critical_pair(i, j);
        }
    }
}

void run_computation() {
    int current_pair_idx = 0;
    while (current_pair_idx < pair_count) {
        CriticalPair pair = critical_pairs[current_pair_idx++];
        Polynomial* p1 = G[pair.p1_idx];
        Polynomial* p2 = G[pair.p2_idx];

        Polynomial* s = calculate_spoly(p1, p2);
        Polynomial* r = reduce_poly(s);
        
        free_poly(s);

        if (r->num_terms > 0) {
            // Normalize to simplify later reductions
            if (r->terms[0].coefficient < -1) r->terms[0].coefficient = -1;
            else if (r->terms[0].coefficient > 1) r->terms[0].coefficient = 1;

            int new_poly_idx = g_size;
            add_poly_to_basis(r);
            for (int i = 0; i < new_poly_idx; i++) {
                add_critical_pair(i, new_poly_idx);
            }
        } else {
            free_poly(r);
        }
    }

    // Accumulate a final result to prevent dead code elimination
    for (int i = 0; i < g_size; i++) {
        for (int j = 0; j < G[i]->num_terms; j++) {
            final_result += G[i]->terms[j].coefficient;
        }
    }
}

void cleanup() {
    for (int i = 0; i < g_size; i++) {
        free_poly(G[i]);
    }
    free(G);
    free(critical_pairs);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%lld\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
