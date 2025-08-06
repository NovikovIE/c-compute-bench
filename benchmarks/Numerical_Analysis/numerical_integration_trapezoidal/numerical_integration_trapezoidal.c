#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- Mersenne Twister (Do Not Modify) ---
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

// Benchmark parameters and data
int g_polynomial_degree;
double g_lower_bound;
double g_upper_bound;
long long g_num_intervals;
double* g_coefficients = NULL;

// Result
double g_integral_result;

// Function to evaluate the polynomial at a given point x using Horner's method
double eval_polynomial(double x) {
    double result = 0.0;
    for (int i = g_polynomial_degree; i >= 0; i--) {
        result = result * x + g_coefficients[i];
    }
    return result;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s polynomial_degree lower_bound upper_bound num_intervals seed\n", argv[0]);
        exit(1);
    }

    g_polynomial_degree = atoi(argv[1]);
    g_lower_bound = atof(argv[2]);
    g_upper_bound = atof(argv[3]);
    g_num_intervals = atoll(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    if (g_polynomial_degree < 0) {
        fprintf(stderr, "Error: Polynomial degree must be non-negative.\n");
        exit(1);
    }
    if (g_num_intervals <= 0) {
        fprintf(stderr, "Error: Number of intervals must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate memory for polynomial coefficients
    g_coefficients = (double*)malloc((g_polynomial_degree + 1) * sizeof(double));
    if (g_coefficients == NULL) {
        fprintf(stderr, "Failed to allocate memory for coefficients\n");
        exit(1);
    }

    // Generate random coefficients between -1.0 and 1.0
    for (int i = 0; i <= g_polynomial_degree; i++) {
        g_coefficients[i] = (double)mt_rand() / (double)UINT32_MAX * 2.0 - 1.0;
    }
}

void run_computation() {
    double h = (g_upper_bound - g_lower_bound) / (double)g_num_intervals;
    double sum = 0.0;

    // The trapezoidal rule sum is h * [ f(x0)/2 + f(x1) + ... + f(xn-1) + f(xn)/2 ]
    // This implementation calculates the sum first, then multiplies by h to improve precision.

    // Add the first and last points, each weighted by 1/2
    sum += eval_polynomial(g_lower_bound);
    sum += eval_polynomial(g_upper_bound);
    sum /= 2.0;

    // Add the interior points, each weighted by 1
    for (long long i = 1; i < g_num_intervals; i++) {
        double x = g_lower_bound + i * h;
        sum += eval_polynomial(x);
    }

    g_integral_result = sum * h;
}

void cleanup() {
    free(g_coefficients);
    g_coefficients = NULL;
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
    printf("%f\n", g_integral_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
