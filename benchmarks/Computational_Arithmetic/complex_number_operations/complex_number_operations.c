#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- START MERSENNE TWISTER (MT19937) ---
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

// ---- BENCHMARK DATA AND PARAMETERS ----

typedef struct {
    double real;
    double imag;
} Complex;

int NUM_OPERATIONS;
Complex *data1;
Complex *data2;

double final_result; // Accumulated result to prevent dead-code elimination

// ---- BENCHMARK FUNCTIONS ----

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_operations> <seed>\n", argv[0]);
        exit(1);
    }

    NUM_OPERATIONS = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if (NUM_OPERATIONS <= 0) {
        fprintf(stderr, "FATAL: num_operations must be a positive integer.\n");
        exit(1);
    }

    mt_seed(seed);

    data1 = (Complex *)malloc(NUM_OPERATIONS * sizeof(Complex));
    data2 = (Complex *)malloc(NUM_OPERATIONS * sizeof(Complex));

    if (data1 == NULL || data2 == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    for (int i = 0; i < NUM_OPERATIONS; i++) {
        data1[i].real = 2.0 * (mt_rand() / 4294967295.0) - 1.0; // Random double in [-1, 1]
        data1[i].imag = 2.0 * (mt_rand() / 4294967295.0) - 1.0;
        data2[i].real = 2.0 * (mt_rand() / 4294967295.0) - 1.0;
        data2[i].imag = 2.0 * (mt_rand() / 4294967295.0) - 1.0;
    }
}

// Helper functions for complex arithmetic. Marked static inline for performance.
static inline Complex complex_multiply(Complex a, Complex b) {
    Complex result;
    result.real = a.real * b.real - a.imag * b.imag;
    result.imag = a.real * b.imag + a.imag * b.real;
    return result;
}

static inline Complex complex_add(Complex a, Complex b) {
    Complex result;
    result.real = a.real + b.real;
    result.imag = a.imag + b.imag;
    return result;
}

void run_computation() {
    Complex accumulator = {1.0, 0.0}; // Start with multiplicative identity

    for (int i = 0; i < NUM_OPERATIONS; i++) {
        // Perform a chain of dependent operations
        accumulator = complex_multiply(accumulator, data1[i]);
        accumulator = complex_add(accumulator, data2[i]);
    }

    // Store the magnitude squared of the final complex number as the result.
    // This is a simple, deterministic value derived from the computation.
    final_result = accumulator.real * accumulator.real + accumulator.imag * accumulator.imag;
}

void cleanup() {
    free(data1);
    free(data2);
}

// ---- MAIN FUNCTION ----

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final result to stdout to ensure computation is not optimized away
    printf("%f\n", final_result);

    // Print the time taken to stderr, as required
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
