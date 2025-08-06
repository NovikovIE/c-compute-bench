#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- BEGIN MERSENNE TWISTER (DO NOT MODIFY) ---
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
// --- END MERSENNE TWISTER ---

// --- BENCHMARK DATA AND PARAMETERS ---
int num_iterations;
double precision;
int poly_degree;

// Each row is a set of coefficients for a polynomial
double **poly_coeffs; 
// Each pair of doubles is a search interval [a, b]
double *search_intervals; 
// Stores the minimum found for each function
double *results;          
// Final accumulated result to prevent dead code elimination
double final_result;

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_iterations> <precision> <poly_degree> <seed>\n", argv[0]);
        exit(1);
    }

    num_iterations = atoi(argv[1]);
    precision = atof(argv[2]);
    poly_degree = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    // Ensure poly_degree is at least 2 and even for robust unimodal functions
    if (poly_degree < 2) {
        poly_degree = 2;
    }
    if (poly_degree % 2 != 0) {
        poly_degree++;
    }
    
    mt_seed(seed);

    poly_coeffs = (double **)malloc(num_iterations * sizeof(double *));
    search_intervals = (double *)malloc(num_iterations * 2 * sizeof(double));
    results = (double *)malloc(num_iterations * sizeof(double));

    if (!poly_coeffs || !search_intervals || !results) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    for (int i = 0; i < num_iterations; ++i) {
        poly_coeffs[i] = (double *)malloc((poly_degree + 1) * sizeof(double));
        if (!poly_coeffs[i]) {
            fprintf(stderr, "FATAL: Memory allocation failed for poly_coeffs[%d].\n", i);
            exit(1);
        }

        // Generate polynomial coefficients
        // We want to create a convex ('u-shaped') function to ensure a single minimum
        // in our search interval. A simple way is to ensure the highest-degree
        // coefficient is positive.
        for (int j = 0; j < poly_degree; ++j) {
            // Random coefficients between -1.0 and 1.0
            poly_coeffs[i][j] = 2.0 * ((double)mt_rand() / (double)UINT32_MAX) - 1.0;
        }
        // Positive coefficient for the highest degree term
        poly_coeffs[i][poly_degree] = (double)mt_rand() / (double)UINT32_MAX;

        // Generate a random search interval [a, b] of fixed width
        double offset = 2.0 * ((double)mt_rand() / (double)UINT32_MAX) - 1.0; // [-1.0, 1.0]
        search_intervals[i * 2] = -2.0 + offset;
        search_intervals[i * 2 + 1] = 2.0 + offset;
    }
}

// Helper: Evaluate a polynomial f(x) using Horner's method for efficiency.
// f(x) = c_d*x^d + c_{d-1}*x^{d-1} + ... + c_1*x + c_0
static double evaluate_poly(double x, const double* coeffs, int degree) {
    double result = 0.0;
    for (int i = degree; i >= 0; --i) {
        result = result * x + coeffs[i];
    }
    return result;
}

// Helper: Find the minimum of a unimodal function on the interval [a, b].
static double find_minimum(double a, double b, const double* coeffs, int degree, double tol) {
    const double inv_phi = (sqrt(5.0) - 1.0) / 2.0; // 1/phi

    double c = b - inv_phi * (b - a);
    double d = a + inv_phi * (b - a);
    double yc = evaluate_poly(c, coeffs, degree);
    double yd = evaluate_poly(d, coeffs, degree);

    while ((b - a) > tol) {
        if (yc < yd) {
            b = d;
            d = c;
            yd = yc;
            c = b - inv_phi * (b - a);
            yc = evaluate_poly(c, coeffs, degree);
        } else {
            a = c;
            c = d;
            yc = yd;
            d = a + inv_phi * (b - a);
            yd = evaluate_poly(d, coeffs, degree);
        }
    }
    return (a + b) / 2.0;
}

void run_computation() {
    for (int i = 0; i < num_iterations; ++i) {
        double a = search_intervals[i * 2];
        double b = search_intervals[i * 2 + 1];
        const double* coeffs = poly_coeffs[i];
        results[i] = find_minimum(a, b, coeffs, poly_degree, precision);
    }

    final_result = 0.0;
    for (int i = 0; i < num_iterations; ++i) {
        final_result += results[i];
    }
}

void cleanup() {
    if (poly_coeffs) {
        for (int i = 0; i < num_iterations; ++i) {
            free(poly_coeffs[i]);
        }
        free(poly_coeffs);
    }
    free(search_intervals);
    free(results);
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
    printf("%f\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
