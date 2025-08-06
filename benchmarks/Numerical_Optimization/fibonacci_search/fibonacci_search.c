#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

//-- BEGIN MERSENNE TWISTER ---------------------------------------------------
// Included verbatim as required
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
//-- END MERSENNE TWISTER -----------------------------------------------------

// Benchmark parameters
static int num_iterations;
static double precision;

// Data for the problems to be solved
typedef struct {
    double a; // Search interval start
    double b; // Search interval end
    double c; // Parameter for the objective function f(x) = (x-c)^2
} Problem;

static Problem* problems;
static double final_result;

// Pre-computed Fibonacci numbers for the search algorithm
#define MAX_FIB_NUM 93 // F_93 is the largest that fits in unsigned long long
static unsigned long long* fib_numbers;

// The function to be minimized. A simple quadratic for demonstration.
static inline double objective_function(double x, double c) {
    double diff = x - c;
    return diff * diff;
}

// Core Fibonacci search implementation
static double fibonacci_search(double a, double b, double c, double tol) {
    int k = 0;
    while (fib_numbers[k] <= (b - a) / tol) {
        k++;
        if (k >= MAX_FIB_NUM) {
            k = MAX_FIB_NUM - 1; // Cap to prevent overflow
            break;
        }
    }

    double x1 = a + ((double)fib_numbers[k - 2] / fib_numbers[k]) * (b - a);
    double x2 = a + ((double)fib_numbers[k - 1] / fib_numbers[k]) * (b - a);
    double f1 = objective_function(x1, c);
    double f2 = objective_function(x2, c);

    for (int i = 0; i < k - 2; ++i) {
        if (f1 < f2) {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = a + (b - x2);
            f1 = objective_function(x1, c);
        } else {
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = b - (x1 - a);
            f2 = objective_function(x2, c);
        }
    }

    // Final interval converges, return its midpoint
    return (a + b) / 2.0;
}

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_iterations> <precision> <seed>\n", argv[0]);
        exit(1);
    }

    num_iterations = atoi(argv[1]);
    precision = atof(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    mt_seed(seed);

    // Pre-compute Fibonacci numbers
    fib_numbers = (unsigned long long*)malloc(MAX_FIB_NUM * sizeof(unsigned long long));
    if (fib_numbers == NULL) {
        fprintf(stderr, "Failed to allocate memory for Fibonacci numbers.\n");
        exit(1);
    }
    fib_numbers[0] = 0;
    fib_numbers[1] = 1;
    for (int i = 2; i < MAX_FIB_NUM; ++i) {
        fib_numbers[i] = fib_numbers[i - 1] + fib_numbers[i - 2];
    }

    // Allocate and generate problem data
    problems = (Problem*)malloc(num_iterations * sizeof(Problem));
    if (problems == NULL) {
        fprintf(stderr, "Failed to allocate memory for problems.\n");
        free(fib_numbers);
        exit(1);
    }

    for (int i = 0; i < num_iterations; ++i) {
        // Generate a random center for the quadratic function's minimum
        double c = (double)mt_rand() / (double)UINT32_MAX * 200.0 - 100.0; // Range [-100, 100]
        // Generate a random search interval width
        double width = (double)mt_rand() / (double)UINT32_MAX * 40.0 + 2.0; // Width [2, 42]
        
        problems[i].c = c;
        problems[i].a = c - width / 2.0;
        problems[i].b = c + width / 2.0;
    }
}

void run_computation() {
    double total_sum_of_minima = 0.0;
    for (int i = 0; i < num_iterations; ++i) {
        Problem p = problems[i];
        double minimum_found = fibonacci_search(p.a, p.b, p.c, precision);
        total_sum_of_minima += minimum_found;
    }
    final_result = total_sum_of_minima;
}

void cleanup() {
    free(problems);
    free(fib_numbers);
}

int main(int argc, char* argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final accumulated result to stdout to prevent dead code elimination
    printf("%.8f\n", final_result);

    // Print timing information to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
