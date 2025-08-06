#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>
#include <stdbool.h>

// --- START: Mersenne Twister (Do Not Modify) ---
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

// --- Data Structures ---
typedef struct {
    long long coefficient;
    int* exponents; // Array of size g_num_variables
} Term;

typedef struct {
    Term* terms;
    int count;
    int capacity;
} Polynomial;

// --- Global Data ---
int g_num_variables;

Polynomial* g_dividend = NULL;
Polynomial* g_divisor = NULL;
Polynomial* g_quotient = NULL;
Polynomial* g_remainder = NULL;

long long g_final_result = 0;

// --- Function Prototypes ---
void setup_benchmark(int argc, char *argv[]);
void run_computation();
void cleanup();

static Polynomial* create_polynomial(int capacity);
static void free_polynomial(Polynomial* p);
static void add_term_copying(Polynomial* p, const Term* t);
static int compare_terms(const void* a, const void* b);
static void sort_polynomial(Polynomial* p);
static void generate_random_polynomial(Polynomial* p, int num_terms, int max_degree);
static Polynomial* multiply_term_by_poly(const Term* t, const Polynomial* p);
static Polynomial* subtract_polys(const Polynomial* p1, const Polynomial* p2);
static void ensure_unique_terms(Polynomial* p);


// --- Main and Benchmark Functions ---

int main(int argc, char *argv[]) {
    struct timespec start, end;

    if (argc != 7) {
        fprintf(stderr, "Usage: %s num_variables dividend_degree divisor_degree dividend_num_terms divisor_num_terms seed\n", argv[0]);
        return 1;
    }

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%lld\n", g_final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}

void setup_benchmark(int argc, char *argv[]) {
    // Parse arguments
    g_num_variables = atoi(argv[1]);
    int dividend_degree = atoi(argv[2]);
    int divisor_degree = atoi(argv[3]);
    int dividend_num_terms = atoi(argv[4]);
    int divisor_num_terms = atoi(argv[5]);
    uint32_t seed = (uint32_t)strtoul(argv[6], NULL, 10);
    
    mt_seed(seed);

    // Create and generate divisor
    g_divisor = create_polynomial(divisor_num_terms);
    generate_random_polynomial(g_divisor, divisor_num_terms, divisor_degree);
    sort_polynomial(g_divisor);
    ensure_unique_terms(g_divisor);
    if (g_divisor->count > 0) {
        // For stable division, make leading term coefficient 1.
        g_divisor->terms[0].coefficient = 1;
    } else {
        fprintf(stderr, "FATAL: Divisor has zero terms after generation.\n");
        exit(1);
    }

    // Create and generate dividend
    g_dividend = create_polynomial(dividend_num_terms);
    generate_random_polynomial(g_dividend, dividend_num_terms, dividend_degree);
    sort_polynomial(g_dividend);
    ensure_unique_terms(g_dividend);
    if (g_dividend->count > 0) {
        // Ensure dividend's leading term is divisible by divisor's leading term.
        Term* lt_div = &g_divisor->terms[0];
        Term* lt_divd = &g_dividend->terms[0];
        for (int i = 0; i < g_num_variables; ++i) {
            lt_divd->exponents[i] = lt_div->exponents[i] + (mt_rand() % 3);
        }
        sort_polynomial(g_dividend);
    } else {
        fprintf(stderr, "FATAL: Dividend has zero terms after generation.\n");
        exit(1);
    }

    // Initialize quotient and remainder polynomials
    g_quotient = create_polynomial(16);
    g_remainder = create_polynomial(g_dividend->count);
    for(int i=0; i < g_dividend->count; ++i) {
        add_term_copying(g_remainder, &g_dividend->terms[i]);
    }
}

void run_computation() {
    while (g_remainder->count > 0) {
        Term* lt_rem = &g_remainder->terms[0];
        Term* lt_div = &g_divisor->terms[0];
        
        bool divisible = true;
        if (g_divisor->count == 0 || lt_rem->coefficient == 0) {
            divisible = false;
        }
        if (divisible) {
            for (int i = 0; i < g_num_variables; ++i) {
                if (lt_rem->exponents[i] < lt_div->exponents[i]) {
                    divisible = false;
                    break;
                }
            }
        }

        if (!divisible) {
            break; // Division process is finished
        }

        Term t; // The next term of the quotient
        t.coefficient = lt_rem->coefficient / lt_div->coefficient;
        t.exponents = (int*)malloc(g_num_variables * sizeof(int));
        for (int i = 0; i < g_num_variables; ++i) {
            t.exponents[i] = lt_rem->exponents[i] - lt_div->exponents[i];
        }
        
        add_term_copying(g_quotient, &t);

        Polynomial* product = multiply_term_by_poly(&t, g_divisor);
        Polynomial* next_remainder = subtract_polys(g_remainder, product);
        
        free(t.exponents);
        free_polynomial(product);
        free_polynomial(g_remainder);
        g_remainder = next_remainder;
    }

    // Accumulate a final result to prevent dead code elimination
    long long sum = 0;
    for (int i = 0; i < g_quotient->count; ++i) {
        sum += g_quotient->terms[i].coefficient;
    }
    for (int i = 0; i < g_remainder->count; ++i) {
        sum += g_remainder->terms[i].coefficient;
    }
    g_final_result = sum;
}

void cleanup() {
    free_polynomial(g_dividend);
    free_polynomial(g_divisor);
    free_polynomial(g_quotient);
    free_polynomial(g_remainder);
}

// --- Helper Function Implementations ---

static Polynomial* create_polynomial(int capacity) {
    Polynomial* p = (Polynomial*)malloc(sizeof(Polynomial));
    if (!p) exit(1);
    p->count = 0;
    p->capacity = capacity;
    if (capacity > 0) {
        p->terms = (Term*)malloc(capacity * sizeof(Term));
        if (!p->terms) exit(1);
    } else {
        p->terms = NULL;
    }
    return p;
}

static void free_polynomial(Polynomial* p) {
    if (p) {
        for (int i = 0; i < p->count; ++i) {
            free(p->terms[i].exponents);
        }
        free(p->terms);
        free(p);
    }
}

static void add_term_copying(Polynomial* p, const Term* t) {
    if (p->count >= p->capacity) {
        p->capacity = p->capacity == 0 ? 16 : p->capacity * 2;
        p->terms = (Term*)realloc(p->terms, p->capacity * sizeof(Term));
        if (!p->terms) exit(1);
    }
    int i = p->count;
    p->terms[i].coefficient = t->coefficient;
    p->terms[i].exponents = (int*)malloc(g_num_variables * sizeof(int));
    if (!p->terms[i].exponents) exit(1);
    memcpy(p->terms[i].exponents, t->exponents, g_num_variables * sizeof(int));
    p->count++;
}

static int compare_terms(const void* a, const void* b) {
    const Term* t1 = (const Term*)a;
    const Term* t2 = (const Term*)b;
    for (int i = 0; i < g_num_variables; ++i) {
        if (t1->exponents[i] > t2->exponents[i]) return -1; // Lexicographical: higher exponent comes first
        if (t1->exponents[i] < t2->exponents[i]) return 1;
    }
    return 0;
}

static void sort_polynomial(Polynomial* p) {
    if (p->count > 1) {
        qsort(p->terms, p->count, sizeof(Term), compare_terms);
    }
}

static void generate_random_polynomial(Polynomial* p, int num_terms, int max_degree) {
    for (int i = 0; i < num_terms; ++i) {
        Term t;
        long long coeff = (mt_rand() % 200) - 99; // Range [-99, 100]
        t.coefficient = (coeff == 0) ? 1 : coeff;
        t.exponents = (int*)malloc(g_num_variables * sizeof(int));
        if (!t.exponents) exit(1);

        int current_degree = 0;
        for (int j = 0; j < g_num_variables - 1; ++j) {
            if (current_degree < max_degree) {
                int exp = mt_rand() % (max_degree - current_degree + 1);
                t.exponents[j] = exp;
                current_degree += exp;
            } else {
                t.exponents[j] = 0;
            }
        }
        t.exponents[g_num_variables - 1] = max_degree - current_degree;
        
        // The add_term_copying will manage memory for us
        add_term_copying(p, &t);
        free(t.exponents); // free the temporary exponents array
    }
}

static void ensure_unique_terms(Polynomial* p) {
    if (p->count <= 1) return;
    int new_count = 1;
    for (int i = 1; i < p->count; i++) {
        if (compare_terms(&p->terms[new_count - 1], &p->terms[i]) == 0) {
            // Same exponents, combine coefficients
            p->terms[new_count-1].coefficient += p->terms[i].coefficient;
            free(p->terms[i].exponents); // Free the duplicate term's exponents
        } else {
            // Different term, move it to the next unique spot
            if(new_count != i) {
                p->terms[new_count] = p->terms[i];
            }
            new_count++;
        }
    }
    // Remove zero-coefficient terms that might have resulted from merging
    int final_count = 0;
    for(int i = 0; i < new_count; ++i) {
        if (p->terms[i].coefficient != 0) {
            if(final_count != i) {
                p->terms[final_count] = p->terms[i];
            }
            final_count++;
        } else {
            free(p->terms[i].exponents);
        }
    }
    p->count = final_count;
}

static Polynomial* multiply_term_by_poly(const Term* t, const Polynomial* p) {
    Polynomial* result = create_polynomial(p->count);
    for (int i = 0; i < p->count; ++i) {
        Term new_term;
        new_term.coefficient = t->coefficient * p->terms[i].coefficient;
        new_term.exponents = (int*)malloc(g_num_variables * sizeof(int));
        if (!new_term.exponents) exit(1);
        for (int j = 0; j < g_num_variables; ++j) {
            new_term.exponents[j] = t->exponents[j] + p->terms[i].exponents[j];
        }
        add_term_copying(result, &new_term);
        free(new_term.exponents);
    }
    return result;
}


static Polynomial* subtract_polys(const Polynomial* p1, const Polynomial* p2) {
    Polynomial* result = create_polynomial(p1->count + p2->count);
    int i = 0, j = 0;

    while (i < p1->count && j < p2->count) {
        int cmp = compare_terms(&p1->terms[i], &p2->terms[j]);
        if (cmp > 0) { // p1 term is larger, copy it
            add_term_copying(result, &p1->terms[i++]);
        } else if (cmp < 0) { // p2 term is larger, copy its negation
            Term t2_neg = p2->terms[j++];
            t2_neg.coefficient = -t2_neg.coefficient;
            add_term_copying(result, &t2_neg);
        } else { // exponents are the same
            long long new_coeff = p1->terms[i].coefficient - p2->terms[j].coefficient;
            if (new_coeff != 0) {
                Term new_term = p1->terms[i];
                new_term.coefficient = new_coeff;
                add_term_copying(result, &new_term);
            }
            i++; j++;
        }
    }

    // Append remaining terms
    while (i < p1->count) {
        add_term_copying(result, &p1->terms[i++]);
    }
    while (j < p2->count) {
        Term t2_neg = p2->terms[j++];
        t2_neg.coefficient = -t2_neg.coefficient;
        add_term_copying(result, &t2_neg);
    }
    return result;
}