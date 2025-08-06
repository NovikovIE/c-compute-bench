#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <string.h>

// --- Mersenne Twister (DO NOT MODIFY) ---
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
// --- End Mersenne Twister ---

// --- Benchmark Data Structures ---
typedef struct {
    long long coefficient;
    int *powers;
} Term;

typedef struct {
    Term *terms;
    size_t count;
    size_t capacity;
} Polynomial;

// --- Global State ---
int g_num_variables;
int g_initial_terms;
int g_expansion_power;

Polynomial *g_base_poly = NULL;
Polynomial *g_result_poly = NULL;
long long g_checksum = 0;

// --- Helper Functions Declaration ---
void free_poly(Polynomial *p);

// --- Core Functions ---

// Comparison function for qsort, lexicographically compares powers
int compare_terms(const void *a, const void *b) {
    Term *term_a = (Term *)a;
    Term *term_b = (Term *)b;
    for (int i = 0; i < g_num_variables; ++i) {
        if (term_a->powers[i] < term_b->powers[i]) return -1;
        if (term_a->powers[i] > term_b->powers[i]) return 1;
    }
    return 0;
}

// Simplifies a polynomial by sorting and merging like terms
void simplify_poly(Polynomial *p) {
    if (p->count <= 1) return;

    qsort(p->terms, p->count, sizeof(Term), compare_terms);

    size_t write_idx = 0;
    for (size_t read_idx = 1; read_idx < p->count; ++read_idx) {
        if (compare_terms(&p->terms[write_idx], &p->terms[read_idx]) == 0) {
            p->terms[write_idx].coefficient += p->terms[read_idx].coefficient;
            free(p->terms[read_idx].powers); // Free the redundant powers array
        } else {
            write_idx++;
            if (write_idx != read_idx) {
                p->terms[write_idx] = p->terms[read_idx];
            }
        }
    }
    p->count = write_idx + 1;
}

// Multiplies two polynomials, returns a new, simplified polynomial
Polynomial* multiply_polys(const Polynomial *p1, const Polynomial *p2) {
    size_t max_terms = p1->count * p2->count;
    if (max_terms == 0) {
        Polynomial* empty_poly = (Polynomial*)calloc(1, sizeof(Polynomial));
        return empty_poly;
    }

    Polynomial* result = (Polynomial*)malloc(sizeof(Polynomial));
    result->capacity = max_terms;
    result->terms = (Term*)malloc(result->capacity * sizeof(Term));
    result->count = 0;

    for (size_t i = 0; i < p1->count; ++i) {
        for (size_t j = 0; j < p2->count; ++j) {
            Term *t1 = &p1->terms[i];
            Term *t2 = &p2->terms[j];
            Term *new_term = &result->terms[result->count++];

            new_term->coefficient = t1->coefficient * t2->coefficient;
            new_term->powers = (int*)malloc(g_num_variables * sizeof(int));
            for (int k = 0; k < g_num_variables; ++k) {
                new_term->powers[k] = t1->powers[k] + t2->powers[k];
            }
        }
    }

    simplify_poly(result);
    return result;
}


// Deep copies a polynomial
Polynomial* deep_copy_poly(const Polynomial *src) {
    Polynomial* dest = (Polynomial*)malloc(sizeof(Polynomial));
    dest->count = src->count;
    dest->capacity = src->count;
    dest->terms = (Term*)malloc(dest->capacity * sizeof(Term));

    for (size_t i = 0; i < src->count; ++i) {
        dest->terms[i].coefficient = src->terms[i].coefficient;
        dest->terms[i].powers = (int*)malloc(g_num_variables * sizeof(int));
        memcpy(dest->terms[i].powers, src->terms[i].powers, g_num_variables * sizeof(int));
    }
    return dest;
}


void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_variables initial_terms expansion_power seed\n", argv[0]);
        exit(1);
    }
    g_num_variables = atoi(argv[1]);
    g_initial_terms = atoi(argv[2]);
    g_expansion_power = atoi(argv[3]);
    uint32_t seed = (uint32_t)strtoul(argv[4], NULL, 10);
    mt_seed(seed);
    
    if (g_num_variables <= 0 || g_initial_terms <= 0 || g_expansion_power < 0) {
        fprintf(stderr, "FATAL: Invalid parameters.\n");
        exit(1);
    }

    g_base_poly = (Polynomial*)malloc(sizeof(Polynomial));
    g_base_poly->count = g_initial_terms;
    g_base_poly->capacity = g_initial_terms;
    g_base_poly->terms = (Term*)malloc(g_initial_terms * sizeof(Term));

    for (int i = 0; i < g_initial_terms; ++i) {
        g_base_poly->terms[i].coefficient = (mt_rand() % 10) + 1;
        g_base_poly->terms[i].powers = (int*)malloc(g_num_variables * sizeof(int));
        for (int j = 0; j < g_num_variables; ++j) {
            g_base_poly->terms[i].powers[j] = mt_rand() % 3;
        }
    }
    simplify_poly(g_base_poly);
}

void run_computation() {
    if (g_expansion_power == 0) {
        g_result_poly = (Polynomial*)calloc(1, sizeof(Polynomial));
        g_result_poly->capacity = 1;
        g_result_poly->count = 1;
        g_result_poly->terms = (Term*)calloc(1, sizeof(Term));
        g_result_poly->terms[0].coefficient = 1;
        g_result_poly->terms[0].powers = (int*)calloc(g_num_variables, sizeof(int));
    } else {
        Polynomial *current_poly = deep_copy_poly(g_base_poly);
        for (int i = 1; i < g_expansion_power; ++i) {
            Polynomial *next_poly = multiply_polys(current_poly, g_base_poly);
            free_poly(current_poly);
            current_poly = next_poly;
        }
        g_result_poly = current_poly;
    }
    
    g_checksum = 0;
    for (size_t i = 0; i < g_result_poly->count; ++i) {
        g_checksum += g_result_poly->terms[i].coefficient;
    }
}

void cleanup() {
    free_poly(g_base_poly);
    g_base_poly = NULL;
    free_poly(g_result_poly);
    g_result_poly = NULL;
}


void free_poly(Polynomial *p) {
    if (p == NULL) return;
    for (size_t i = 0; i < p->count; ++i) {
        free(p->terms[i].powers);
    }
    free(p->terms);
    free(p);
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
    printf("%lld\n", g_checksum);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
