#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

// --- Mersenne Twister (MT19937) --- (DO NOT MODIFY)
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
// --- End of Mersenne Twister ---

// --- Benchmark Specifics ---

// Parameters
static int num_variables;
static int polynomial_degree;
static int max_coefficient_bits;
static int num_candidate_factors;

// Data Structures
typedef struct {
    int64_t coefficient;
    int* exponents; // Array of size num_variables
} Term;

typedef struct {
    Term* terms;
    size_t num_terms;
} Polynomial;

// Global Data
static Polynomial* g_candidate_factors = NULL;
static Polynomial* g_final_poly = NULL;
static long long g_final_result = 0;

// --- Helper Functions ---

int64_t rand_64() {
    uint64_t r1 = mt_rand();
    uint64_t r2 = mt_rand();
    return (int64_t)((r1 << 32) | r2);
}

int compare_exponents(const int* exp1, const int* exp2) {
    return memcmp(exp1, exp2, num_variables * sizeof(int));
}

int compare_terms(const void* a, const void* b) {
    const Term* term_a = (const Term*)a;
    const Term* term_b = (const Term*)b;
    return compare_exponents(term_a->exponents, term_b->exponents);
}

void free_polynomial_content(Polynomial* p) {
    if (!p) return;
    for (size_t i = 0; i < p->num_terms; ++i) {
        free(p->terms[i].exponents);
    }
    free(p->terms);
}

Polynomial* deep_copy_polynomial(const Polynomial* src) {
    Polynomial* dst = (Polynomial*)malloc(sizeof(Polynomial));
    if (!dst) { perror("malloc"); exit(1); }
    dst->num_terms = src->num_terms;
    dst->terms = (Term*)malloc(dst->num_terms * sizeof(Term));
    if (!dst->terms) { perror("malloc"); exit(1); }

    for (size_t i = 0; i < src->num_terms; ++i) {
        dst->terms[i].coefficient = src->terms[i].coefficient;
        dst->terms[i].exponents = (int*)malloc(num_variables * sizeof(int));
        if (!dst->terms[i].exponents) { perror("malloc"); exit(1); }
        memcpy(dst->terms[i].exponents, src->terms[i].exponents, num_variables * sizeof(int));
    }
    return dst;
}

Polynomial* polynomial_multiply(const Polynomial* p1, const Polynomial* p2) {
    if (p1->num_terms == 0 || p2->num_terms == 0) {
        Polynomial* result = (Polynomial*)calloc(1, sizeof(Polynomial));
        if (!result) { perror("calloc"); exit(1); }
        return result;
    }

    size_t unsimplified_size = p1->num_terms * p2->num_terms;
    Term* temp_terms = (Term*)malloc(unsimplified_size * sizeof(Term));
    if (!temp_terms) { perror("malloc"); exit(1); }

    size_t k = 0;
    for (size_t i = 0; i < p1->num_terms; ++i) {
        for (size_t j = 0; j < p2->num_terms; ++j) {
            temp_terms[k].coefficient = p1->terms[i].coefficient * p2->terms[j].coefficient;
            temp_terms[k].exponents = (int*)malloc(num_variables * sizeof(int));
            if (!temp_terms[k].exponents) { perror("malloc"); exit(1); }
            for (int v = 0; v < num_variables; ++v) {
                temp_terms[k].exponents[v] = p1->terms[i].exponents[v] + p2->terms[j].exponents[v];
            }
            k++;
        }
    }

    qsort(temp_terms, unsimplified_size, sizeof(Term), compare_terms);

    size_t simplified_size = 0;
    if (unsimplified_size > 0) {
        simplified_size = 1;
        for (size_t i = 1; i < unsimplified_size; i++) {
            if (compare_exponents(temp_terms[i - 1].exponents, temp_terms[i].exponents) != 0) {
                simplified_size++;
            }
        }
    }

    Polynomial* result = (Polynomial*)malloc(sizeof(Polynomial));
    if (!result) { perror("malloc"); exit(1); }
    result->terms = (Term*)malloc(simplified_size * sizeof(Term));
    if (!result->terms) { perror("malloc"); exit(1); }
    result->num_terms = simplified_size;

    if (unsimplified_size > 0) {
        result->terms[0] = temp_terms[0]; 
        size_t result_idx = 0;
        for (size_t i = 1; i < unsimplified_size; ++i) {
            if (compare_exponents(result->terms[result_idx].exponents, temp_terms[i].exponents) == 0) {
                result->terms[result_idx].coefficient += temp_terms[i].coefficient;
                free(temp_terms[i].exponents);
            } else {
                result_idx++;
                result->terms[result_idx] = temp_terms[i];
            }
        }
    }

    free(temp_terms);
    return result;
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s num_variables polynomial_degree max_coefficient_bits num_candidate_factors seed\n", argv[0]);
        exit(1);
    }
    num_variables = atoi(argv[1]);
    polynomial_degree = atoi(argv[2]);
    max_coefficient_bits = atoi(argv[3]);
    num_candidate_factors = atoi(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    mt_seed(seed);

    g_candidate_factors = (Polynomial*)malloc(num_candidate_factors * sizeof(Polynomial));
    if (!g_candidate_factors) { perror("malloc"); exit(1); }

    int64_t coeff_mask = (max_coefficient_bits >= 64) ? -1LL : (1LL << (max_coefficient_bits - 1)) - 1;
    const int initial_terms = 2;

    for (int i = 0; i < num_candidate_factors; i++) {
        g_candidate_factors[i].num_terms = initial_terms;
        g_candidate_factors[i].terms = (Term*)malloc(initial_terms * sizeof(Term));
        if (!g_candidate_factors[i].terms) { perror("malloc"); exit(1); }

        for (int j = 0; j < initial_terms; j++) {
            Term* term = &g_candidate_factors[i].terms[j];
            term->coefficient = (rand_64() & coeff_mask) + 1;
            if (mt_rand() % 2) term->coefficient *= -1;

            term->exponents = (int*)calloc(num_variables, sizeof(int));
            if (!term->exponents) { perror("calloc"); exit(1); }
            
            int remaining_degree = polynomial_degree;
            for (int v = 0; v < num_variables - 1; ++v) {
                if (remaining_degree > 0) {
                    int exp = mt_rand() % (remaining_degree + 1);
                    term->exponents[v] = exp;
                    remaining_degree -= exp;
                }
            }
            term->exponents[num_variables - 1] = remaining_degree;
        }
    }
}

void run_computation() {
    if (num_candidate_factors == 0) {
        g_final_poly = (Polynomial*)calloc(1, sizeof(Polynomial));
        return;
    }

    g_final_poly = deep_copy_polynomial(&g_candidate_factors[0]);

    for (int i = 1; i < num_candidate_factors; i++) {
        Polynomial* product = polynomial_multiply(g_final_poly, &g_candidate_factors[i]);
        free_polynomial_content(g_final_poly);
        free(g_final_poly);
        g_final_poly = product;
    }

    long long sum = 0;
    for (size_t i = 0; i < g_final_poly->num_terms; i++) {
        sum ^= g_final_poly->terms[i].coefficient;
    }
    g_final_result = sum;
}

void cleanup() {
    if (g_candidate_factors) {
        for (int i = 0; i < num_candidate_factors; ++i) {
            free_polynomial_content(&g_candidate_factors[i]);
        }
        free(g_candidate_factors);
    }

    if (g_final_poly) {
        free_polynomial_content(g_final_poly);
        free(g_final_poly);
    }
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

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
