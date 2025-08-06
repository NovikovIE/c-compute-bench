#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- Mersenne Twister (MT19937) Generator ---
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

// --- Benchmark Globals ---
int num_digits;
int* digits_output;     // To store the computed digits of e
int* spigot_coeffs;     // Internal state for the spigot algorithm
int final_result;       // Accumulated result to prevent dead-code elimination

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_digits> <seed>\n", argv[0]);
        exit(1);
    }

    num_digits = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if (num_digits <= 0) {
        fprintf(stderr, "ERROR: num_digits must be a positive integer.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate memory for the output digits
    digits_output = (int*)malloc(num_digits * sizeof(int));
    if (digits_output == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed for digits_output.\n");
        exit(1);
    }

    // Allocate memory for spigot algorithm coefficients. Size is num_digits + 1.
    spigot_coeffs = (int*)malloc((num_digits + 1) * sizeof(int));
    if (spigot_coeffs == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed for spigot_coeffs.\n");
        free(digits_output);
        exit(1);
    }
}

void run_computation() {
    // This benchmark computes the first 'num_digits' of 'e' using a spigot algorithm.
    // 'e' is approximately 2.71828...
    // The algorithm works by maintaining an array of coefficients in a mixed radix representation.

    const int array_len = num_digits + 1;

    // Initialize all coefficients to 1
    for (int i = 0; i < array_len; i++) {
        spigot_coeffs[i] = 1;
    }

    // The first digit of 'e' is always 2
    digits_output[0] = 2;

    // Compute the remaining num_digits - 1 fractional digits
    for (int k = 1; k < num_digits; k++) {
        int carry = 0;
        // Process coefficients from right to left
        for (int i = num_digits; i >= 1; i--) {
            int temp = spigot_coeffs[i] * 10 + carry;
            // Divisor is (i+1) corresponding to the terms in the series for e-2
            spigot_coeffs[i] = temp % (i + 1);
            carry = temp / (i + 1);
        }
        digits_output[k] = carry;
    }

    // To ensure the computation is not optimized away, we sum all computed digits.
    int sum = 0;
    for (int i = 0; i < num_digits; i++) {
        sum += digits_output[i];
    }
    final_result = sum;
}

void cleanup() {
    free(digits_output);
    free(spigot_coeffs);
}

// --- Main Function ---

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout
    printf("%d\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
