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

// Benchmark-specific data structures and globals
typedef struct {
    double lower;
    double upper;
} Interval;

int num_operations;
Interval* operands_a = NULL;
Interval* operands_b = NULL;
Interval* results = NULL;
double final_result_accumulator = 0.0;

// Helper to generate a random double between 0.0 and 1.0
double rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

// Helper functions for interval multiplication
static inline double min4(double a, double b, double c, double d) {
    return fmin(fmin(a, b), fmin(c, d));
}

static inline double max4(double a, double b, double c, double d) {
    return fmax(fmax(a, b), fmax(c, d));
}

// Interval arithmetic operations
static inline Interval interval_add(Interval x, Interval y) {
    Interval result = {x.lower + y.lower, x.upper + y.upper};
    return result;
}

static inline Interval interval_sub(Interval x, Interval y) {
    Interval result = {x.lower - y.upper, x.upper - y.lower};
    return result;
}

static inline Interval interval_mul(Interval x, Interval y) {
    double p1 = x.lower * y.lower;
    double p2 = x.lower * y.upper;
    double p3 = x.upper * y.lower;
    double p4 = x.upper * y.upper;
    Interval result = {min4(p1, p2, p3, p4), max4(p1, p2, p3, p4)};
    return result;
}

// Assumes y does not contain 0. Our setup guarantees y is strictly positive.
static inline Interval interval_div(Interval x, Interval y) {
    Interval recip_y = {1.0 / y.upper, 1.0 / y.lower};
    return interval_mul(x, recip_y);
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_operations> <seed>\n", argv[0]);
        exit(1);
    }
    
    num_operations = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);
    
    if (num_operations <= 0) {
        fprintf(stderr, "FATAL: num_operations must be > 0\n");
        exit(1);
    }

    mt_seed(seed);

    operands_a = (Interval*)malloc(num_operations * sizeof(Interval));
    operands_b = (Interval*)malloc(num_operations * sizeof(Interval));
    results = (Interval*)malloc(num_operations * sizeof(Interval));

    if (!operands_a || !operands_b || !results) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    for (int i = 0; i < num_operations; i++) {
        // Operand A can be any interval
        double r1 = rand_double() * 200.0 - 100.0; // [-100, 100]
        double r2 = rand_double() * 200.0 - 100.0;
        operands_a[i].lower = fmin(r1, r2);
        operands_a[i].upper = fmax(r1, r2);

        // Operand B must not contain zero for division. We make it strictly positive.
        double r3 = rand_double() * 99.0 + 1.0; // [1, 100]
        double r4 = rand_double() * 99.0 + 1.0;
        operands_b[i].lower = fmin(r3, r4);
        operands_b[i].upper = fmax(r3, r4);
    }
}

void run_computation() {
    double accumulator = 0.0;
    for (int i = 0; i < num_operations; i++) {
        Interval a = operands_a[i];
        Interval b = operands_b[i];
        
        // Perform a chain of operations
        Interval res_add = interval_add(a, b);
        Interval res_sub = interval_sub(a, b);
        Interval res_mul = interval_mul(res_add, res_sub);
        Interval res_div = interval_div(res_mul, b);

        results[i] = res_div;
        accumulator += res_div.lower + res_div.upper;
    }
    final_result_accumulator = accumulator;
}

void cleanup() {
    free(operands_a);
    free(operands_b);
    free(results);
    operands_a = NULL;
    operands_b = NULL;
    results = NULL;
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout to prevent dead code elimination
    printf("%f\n", final_result_accumulator);

    // Print timing information to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
