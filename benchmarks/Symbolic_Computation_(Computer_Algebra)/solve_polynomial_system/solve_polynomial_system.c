#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

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

// --- BENCHMARK DATA STRUCTURES ---

typedef struct {
    long long coefficient;
    int* exponents; // Array of size NUM_VARIABLES
} Term;

typedef struct {
    Term* terms;
    int num_terms;
} Polynomial;

// --- GLOBAL VARIABLES ---

// Parameters
int NUM_EQUATIONS;
int NUM_VARIABLES;
int MAX_DEGREE;

// Data
Polynomial** system_of_equations = NULL;
Polynomial* final_result_polynomial = NULL;
long long final_checksum = 0;

// --- FORWARD DECLARATIONS ---
void free_polynomial(Polynomial* p);
int term_comparator(const void* a, const void* b);
void simplify_polynomial(Polynomial* p);

// --- BENCHMARK IMPLEMENTATION ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_equations> <num_variables> <max_degree> <seed>\n", argv[0]);
        exit(1);
    }

    NUM_EQUATIONS = atoi(argv[1]);
    NUM_VARIABLES = atoi(argv[2]);
    MAX_DEGREE = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    mt_seed(seed);

    system_of_equations = (Polynomial**)malloc(NUM_EQUATIONS * sizeof(Polynomial*));
    if (!system_of_equations) {
        perror("Failed to allocate system of equations");
        exit(1);
    }

    for (int i = 0; i < NUM_EQUATIONS; i++) {
        int num_terms = 2 + (mt_rand() % 2); // 2 to 3 terms
        Polynomial* p = (Polynomial*)malloc(sizeof(Polynomial));
        if (!p) { perror("malloc failed"); exit(1); }
        
        p->num_terms = num_terms;
        p->terms = (Term*)malloc(num_terms * sizeof(Term));
        if (!p->terms) { perror("malloc failed"); exit(1); }

        for (int j = 0; j < num_terms; j++) {
            p->terms[j].coefficient = (mt_rand() % 10) - 4; // Coeffs from -4 to 5
            if(p->terms[j].coefficient == 0) p->terms[j].coefficient = 1;

            p->terms[j].exponents = (int*)malloc(NUM_VARIABLES * sizeof(int));
            if (!p->terms[j].exponents) { perror("malloc failed"); exit(1); }

            int current_degree = 0;
            for (int k = 0; k < NUM_VARIABLES; k++) {
                if (current_degree < MAX_DEGREE) {
                    int exponent = mt_rand() % (MAX_DEGREE - current_degree + 1);
                    p->terms[j].exponents[k] = exponent;
                    current_degree += exponent;
                } else {
                    p->terms[j].exponents[k] = 0;
                }
            }
        }
        system_of_equations[i] = p;
    }
}

Polynomial* copy_polynomial(const Polynomial* p_in) {
    Polynomial* p_out = (Polynomial*)malloc(sizeof(Polynomial));
    if (!p_out) { perror("malloc failed"); exit(1); }

    p_out->num_terms = p_in->num_terms;
    p_out->terms = (Term*)malloc(p_out->num_terms * sizeof(Term));
    if (!p_out->terms) { perror("malloc failed"); exit(1); }

    for (int i = 0; i < p_out->num_terms; i++) {
        p_out->terms[i].coefficient = p_in->terms[i].coefficient;
        p_out->terms[i].exponents = (int*)malloc(NUM_VARIABLES * sizeof(int));
        if (!p_out->terms[i].exponents) { perror("malloc failed"); exit(1); }
        for (int j = 0; j < NUM_VARIABLES; j++) {
            p_out->terms[i].exponents[j] = p_in->terms[i].exponents[j];
        }
    }
    return p_out;
}

Polynomial* multiply_polynomials(const Polynomial* p1, const Polynomial* p2) {
    if (p1->num_terms == 0 || p2->num_terms == 0) {
         Polynomial* p_res = (Polynomial*)malloc(sizeof(Polynomial));
         if (!p_res) { perror("malloc failed"); exit(1); }
         p_res->num_terms = 0;
         p_res->terms = NULL;
         return p_res;
    }
    long long max_terms = (long long)p1->num_terms * p2->num_terms;
    Polynomial* p_res = (Polynomial*)malloc(sizeof(Polynomial));
    if (!p_res) { perror("malloc failed"); exit(1); }
    
    p_res->num_terms = max_terms;
    p_res->terms = (Term*)malloc(max_terms * sizeof(Term));
    if (!p_res->terms) { perror("malloc failed"); exit(1); }

    int current_term = 0;
    for (int i = 0; i < p1->num_terms; i++) {
        for (int j = 0; j < p2->num_terms; j++) {
            p_res->terms[current_term].coefficient = p1->terms[i].coefficient * p2->terms[j].coefficient;
            p_res->terms[current_term].exponents = (int*)malloc(NUM_VARIABLES * sizeof(int));
            if (!p_res->terms[current_term].exponents) { perror("malloc failed"); exit(1); }

            for (int k = 0; k < NUM_VARIABLES; k++) {
                p_res->terms[current_term].exponents[k] = p1->terms[i].exponents[k] + p2->terms[j].exponents[k];
            }
            current_term++;
        }
    }
    p_res->num_terms = current_term;
    return p_res;
}

void run_computation() {
    if (NUM_EQUATIONS == 0) {
        final_checksum = 0;
        return;
    }

    Polynomial* result = copy_polynomial(system_of_equations[0]);

    for (int i = 1; i < NUM_EQUATIONS; i++) {
        Polynomial* product = multiply_polynomials(result, system_of_equations[i]);
        free_polynomial(result);
        simplify_polynomial(product);
        result = product;
    }

    final_result_polynomial = result;
    for (int i = 0; i < final_result_polynomial->num_terms; i++) {
        final_checksum += final_result_polynomial->terms[i].coefficient;
    }
}

void cleanup() {
    if (system_of_equations) {
        for (int i = 0; i < NUM_EQUATIONS; i++) {
            free_polynomial(system_of_equations[i]);
        }
        free(system_of_equations);
    }
    free_polynomial(final_result_polynomial);
}

// --- HELPER & UTILITY FUNCTIONS ---

void free_polynomial(Polynomial* p) {
    if (!p) return;
    if (p->terms) {
        for (int i = 0; i < p->num_terms; i++) {
            if(p->terms[i].exponents) free(p->terms[i].exponents);
        }
        free(p->terms);
    }
    free(p);
}

int term_comparator(const void* a, const void* b) {
    Term* termA = (Term*)a;
    Term* termB = (Term*)b;
    for (int i = 0; i < NUM_VARIABLES; ++i) {
        if (termA->exponents[i] < termB->exponents[i]) return 1;
        if (termA->exponents[i] > termB->exponents[i]) return -1;
    }
    return 0;
}

void simplify_polynomial(Polynomial* p) {
    if (p->num_terms <= 1) return;

    qsort(p->terms, p->num_terms, sizeof(Term), term_comparator);

    int write_idx = 0;
    for (int read_idx = 1; read_idx < p->num_terms; read_idx++) {
        if (term_comparator(&p->terms[read_idx], &p->terms[write_idx]) == 0) {
            p->terms[write_idx].coefficient += p->terms[read_idx].coefficient;
            free(p->terms[read_idx].exponents);
            p->terms[read_idx].exponents = NULL;
        } else {
            write_idx++;
            if (write_idx != read_idx) {
                p->terms[write_idx] = p->terms[read_idx];
            }
        }
    }
    int combined_count = write_idx + 1;

    int final_count = 0;
    for (int i = 0; i < combined_count; i++) {
        if (p->terms[i].coefficient != 0) {
            if (final_count != i) {
                p->terms[final_count] = p->terms[i];
            }
            final_count++;
        } else {
            free(p->terms[i].exponents);
        }
    }
    
    p->num_terms = final_count;
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

    printf("%lld\n", final_checksum);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
