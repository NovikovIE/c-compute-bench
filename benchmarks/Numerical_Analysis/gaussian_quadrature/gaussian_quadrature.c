#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <string.h>
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

// Define a function pointer type for the integrand
typedef double (*IntegrandFunc)(double);

// Global variables for benchmark data
static IntegrandFunc integrand_func;
static double lower_bound;
static double upper_bound;
static int num_points;
static double* weights;
static double* abscissas;

// Global variable for the final result
static double final_result;

// Integrand function implementations
double func_poly(double x) {
    return x * x * x - 2 * x * x + 5 * x - 10;
}

double func_sin(double x) {
    return sin(x);
}

double func_exp(double x) {
    return exp(x);
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <function_str> <lower_bound> <upper_bound> <num_points> <seed>\n", argv[0]);
        exit(1);
    }

    char* function_str = argv[1];
    lower_bound = strtod(argv[2], NULL);
    upper_bound = strtod(argv[3], NULL);
    num_points = atoi(argv[4]);
    uint32_t seed = (uint32_t)strtoul(argv[5], NULL, 10);

    if (strcmp(function_str, "poly") == 0) {
        integrand_func = func_poly;
    } else if (strcmp(function_str, "sin") == 0) {
        integrand_func = func_sin;
    } else if (strcmp(function_str, "exp") == 0) {
        integrand_func = func_exp;
    } else {
        fprintf(stderr, "Error: Invalid function_str '%s'. Choose from 'poly', 'sin', 'exp'.\n", function_str);
        exit(1);
    }

    if (num_points <= 0) {
        fprintf(stderr, "Error: num_points must be a positive integer.\n");
        exit(1);
    }

    mt_seed(seed);

    weights = (double*)malloc(num_points * sizeof(double));
    abscissas = (double*)malloc(num_points * sizeof(double));

    if (!weights || !abscissas) {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        exit(1);
    }

    // Generate synthetic (but deterministic) weights and abscissas for the [-1, 1] interval.
    // These are based on Chebyshev nodes, not true Gauss-Legendre points, but serve to
    // create a realistic computational load.
    for (int i = 0; i < num_points; ++i) {
        abscissas[i] = cos((2.0 * i + 1.0) * M_PI / (2.0 * num_points));
        // A simplified weight; real weights are non-uniform. For benchmark purposes,
        // this is sufficient to stress the FPU.
        weights[i] = M_PI / num_points;
    }
}

void run_computation() {
    double sum = 0.0;
    double transform_mult = (upper_bound - lower_bound) / 2.0;
    double transform_add = (upper_bound + lower_bound) / 2.0;

    for (int i = 0; i < num_points; ++i) {
        // Map the abscissa from [-1, 1] to [lower_bound, upper_bound]
        double x = transform_mult * abscissas[i] + transform_add;
        sum += weights[i] * integrand_func(x);
    }

    final_result = transform_mult * sum;
}

void cleanup() {
    if (weights) free(weights);
    if (abscissas) free(abscissas);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final result to stdout
    printf("%.10f\n", final_result);

    // Print the execution time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
