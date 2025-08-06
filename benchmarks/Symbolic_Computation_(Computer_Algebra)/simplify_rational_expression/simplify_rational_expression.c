#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>
#include <stdbool.h>

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

// --- BENCHMARK DATA STRUCTURES ---

// Represents a single term in a polynomial, e.g., 5 * x^2 * y^3
typedef struct {
    int coefficient;  // e.g., 5
    int* exponents;   // Array of exponents for each variable, e.g., [2, 3]
} Term;

// Represents a multivariate polynomial as a sum of terms
typedef struct {
    Term* terms;
    int num_terms;
} Polynomial;

// Holds all benchmark parameters and data pointers
typedef struct {
    int num_variables;
    int numerator_degree;
    int denominator_degree;
    int numerator_terms;
    int denominator_terms;
    uint32_t seed;

    Polynomial numerator;
    Polynomial denominator;
} BenchmarkData;

// --- GLOBAL VARIABLES ---

// Global pointer to benchmark data, shared by setup, compute, and cleanup
BenchmarkData* g_data = NULL;

// Global variable to store the final result, ensuring it survives cleanup()
long long g_final_result = 0;

// --- HELPER FUNCTIONS ---

// Helper to generate a single polynomial
void generate_polynomial(Polynomial* poly, int num_terms, int degree, int num_variables) {
    poly->num_terms = num_terms;
    poly->terms = (Term*)malloc(num_terms * sizeof(Term));
    if (!poly->terms) {
        perror("malloc failed for terms");
        exit(1);
    }

    for (int i = 0; i < num_terms; ++i) {
        poly->terms[i].coefficient = (mt_rand() % 201) - 100; // Coeff in [-100, 100]
        if (poly->terms[i].coefficient == 0) poly->terms[i].coefficient = 1; // Avoid zero coefficients

        poly->terms[i].exponents = (int*)malloc(num_variables * sizeof(int));
        if (!poly->terms[i].exponents) {
            perror("malloc failed for exponents");
            exit(1);
        }

        // Distribute the total degree among the variables
        int remaining_degree = degree;
        for (int j = 0; j < num_variables - 1; ++j) {
            if (remaining_degree > 0) {
                int exp = mt_rand() % (remaining_degree + 1);
                poly->terms[i].exponents[j] = exp;
                remaining_degree -= exp;
            } else {
                poly->terms[i].exponents[j] = 0;
            }
        }
        poly->terms[i].exponents[num_variables - 1] = remaining_degree;
    }
}

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 7) {
        fprintf(stderr, "Usage: %s num_variables num_degree den_degree num_terms den_terms seed\n", argv[0]);
        exit(1);
    }

    g_data = (BenchmarkData*)malloc(sizeof(BenchmarkData));
    if (!g_data) {
        perror("malloc failed for BenchmarkData");
        exit(1);
    }

    // Parse command-line arguments
    g_data->num_variables = atoi(argv[1]);
    g_data->numerator_degree = atoi(argv[2]);
    g_data->denominator_degree = atoi(argv[3]);
    g_data->numerator_terms = atoi(argv[4]);
    g_data->denominator_terms = atoi(argv[5]);
    g_data->seed = (uint32_t)strtoul(argv[6], NULL, 10);
 
    mt_seed(g_data->seed);

    // Generate the numerator and denominator polynomials
    generate_polynomial(&g_data->numerator, g_data->numerator_terms, g_data->numerator_degree, g_data->num_variables);
    generate_polynomial(&g_data->denominator, g_data->denominator_terms, g_data->denominator_degree, g_data->num_variables);
}

void run_computation() {
    long long cancellation_sum = 0;
    const int num_vars = g_data->num_variables;

    // To avoid matching a denominator term multiple times
    bool* matched_denom_terms = (bool*)calloc(g_data->denominator.num_terms, sizeof(bool));
    if (!matched_denom_terms) {
        perror("calloc failed for matched_denom_terms");
        exit(1);
    }

    // The core computation: find common factors (terms with identical exponents)
    // and sum their coefficients as a proxy for cancellation.
    // This is O(num_terms_N * num_terms_D * num_variables).
    for (int i = 0; i < g_data->numerator.num_terms; ++i) {
        for (int j = 0; j < g_data->denominator.num_terms; ++j) {
            if (matched_denom_terms[j]) {
                continue; // This denominator term is already 'cancelled'
            }

            // Compare exponents to find a common factor
            if (memcmp(g_data->numerator.terms[i].exponents, g_data->denominator.terms[j].exponents, num_vars * sizeof(int)) == 0) {
                 // Found a match. Accumulate and mark as matched.
                 cancellation_sum += g_data->numerator.terms[i].coefficient;
                 cancellation_sum -= g_data->denominator.terms[j].coefficient; // Subtract to add variety
                 matched_denom_terms[j] = true;
                 break; // Move to the next numerator term
            }
        }
    }

    g_final_result = cancellation_sum;
    free(matched_denom_terms);
}

void cleanup() {
    if (!g_data) return;

    // Free numerator data
    for (int i = 0; i < g_data->numerator.num_terms; ++i) {
        free(g_data->numerator.terms[i].exponents);
    }
    free(g_data->numerator.terms);

    // Free denominator data
    for (int i = 0; i < g_data->denominator.num_terms; ++i) {
        free(g_data->denominator.terms[i].exponents);
    }
    free(g_data->denominator.terms);

    // Free the main data structure
    free(g_data);
    g_data = NULL;
}


int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);
    
    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // The result is stored in `g_final_result` to survive cleanup.
    cleanup();

    // Print final result to stdout
    printf("%lld\n", g_final_result);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
