#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

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

// Benchmark-specific constants
#define NUM_RUNS 300000
#define MAX_ITERATIONS 100
#define DERIVATIVE_EPSILON 1e-12

// Global data structure
typedef struct {
    int poly_degree;
    double tolerance;
    double* poly_coeffs;      // Coefficients of the polynomial f(x)
    double* deriv_coeffs;     // Coefficients of the derivative f'(x)
    double* initial_guesses;  // Starting points for each Newton-Raphson run
    double final_result_accumulator;
} BenchmarkData;

static BenchmarkData g_data;

// Evaluates a polynomial p(x) using Horner's method for numerical stability and efficiency.
static double evaluate_poly(const double* coeffs, int degree, double x) {
    double result = 0.0;
    for (int i = degree; i >= 0; --i) {
        result = result * x + coeffs[i];
    }
    return result;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <function_str> <poly_degree> <initial_guess> <tolerance> <seed>\n", argv[0]);
        exit(1);
    }

    // function_str (argv[1]) is unused, as we generate a random polynomial.
    g_data.poly_degree = atoi(argv[2]);
    double first_initial_guess = atof(argv[3]);
    g_data.tolerance = atof(argv[4]);
    uint32_t seed = (uint32_t)strtoul(argv[5], NULL, 10);

    if (g_data.poly_degree <= 0) {
        fprintf(stderr, "FATAL: Polynomial degree must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    g_data.poly_coeffs = (double*)malloc((g_data.poly_degree + 1) * sizeof(double));
    g_data.deriv_coeffs = (double*)malloc(g_data.poly_degree * sizeof(double));
    g_data.initial_guesses = (double*)malloc(NUM_RUNS * sizeof(double));

    if (!g_data.poly_coeffs || !g_data.deriv_coeffs || !g_data.initial_guesses) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Generate polynomial coefficients in the range [-5.0, 5.0]
    for (int i = 0; i <= g_data.poly_degree; ++i) {
        g_data.poly_coeffs[i] = (mt_rand() / (double)UINT32_MAX) * 10.0 - 5.0;
    }
    // Ensure the highest degree coefficient is non-zero to a degree `poly_degree` polynomial.
    if (fabs(g_data.poly_coeffs[g_data.poly_degree]) < 1e-6) {
        g_data.poly_coeffs[g_data.poly_degree] = 1.0;
    }

    // Calculate derivative coefficients
    for (int i = 0; i < g_data.poly_degree; ++i) {
        g_data.deriv_coeffs[i] = g_data.poly_coeffs[i + 1] * (double)(i + 1);
    }

    // Populate initial guesses for each independent run
    g_data.initial_guesses[0] = first_initial_guess; // Use the provided guess for the first run
    for (int i = 1; i < NUM_RUNS; ++i) {
        // Generate other initial guesses in the range [-10.0, 10.0]
        g_data.initial_guesses[i] = (mt_rand() / (double)UINT32_MAX) * 20.0 - 10.0;
    }

    g_data.final_result_accumulator = 0.0;
}

void run_computation() {
    double accumulator = 0.0;
    for (int i = 0; i < NUM_RUNS; ++i) {
        double x = g_data.initial_guesses[i];

        for (int iter = 0; iter < MAX_ITERATIONS; ++iter) {
            double fx = evaluate_poly(g_data.poly_coeffs, g_data.poly_degree, x);

            if (fabs(fx) < g_data.tolerance) {
                break; // Converged
            }

            double fpx = evaluate_poly(g_data.deriv_coeffs, g_data.poly_degree - 1, x);

            if (fabs(fpx) < DERIVATIVE_EPSILON) {
                // Derivative is too small, cannot continue division.
                // Break to avoid NaNs which would spoil the accumulator.
                break;
            }

            x = x - fx / fpx;
        }
        accumulator += x; // Accumulate the found root to prevent dead code elimination
    }
    g_data.final_result_accumulator = accumulator;
}

void cleanup() {
    free(g_data.poly_coeffs);
    free(g_data.deriv_coeffs);
    free(g_data.initial_guesses);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final accumulated result to stdout
    printf("%.10f\n", g_data.final_result_accumulator);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
