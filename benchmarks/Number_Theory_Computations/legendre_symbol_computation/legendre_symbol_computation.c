#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

//
// Program: legendre_symbol_computation
// Theme:   Number Theory Computations
//

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
// --- end Mersenne Twister ---


// --- Benchmark Globals ---
int num_pairs;
long long modulus; // This must be an odd prime for Legendre symbol
long long *a_values;
int *l_symbols;
long long final_result;


// --- Helper Functions ---

// Computes (base^exp) % p using modular exponentiation.
// It uses __int128_t for intermediate products to prevent overflow
// when p is large. This is a GCC/Clang extension.
long long power(long long base, long long exp, long long p) {
    long long res = 1;
    base %= p;
    while (exp > 0) {
        if (exp % 2 == 1) {
            __int128_t temp = (__int128_t)res * base;
            res = (long long)(temp % p);
        }
        __int128_t temp = (__int128_t)base * base;
        base = (long long)(temp % p);
        exp /= 2;
    }
    return res;
}

// Computes the Legendre symbol (a/p) using Euler's criterion.
// p must be an odd prime.
// Returns 1 if a is a quadratic residue mod p.
// Returns -1 if a is a quadratic non-residue mod p.
// Returns 0 if a is divisible by p.
int legendre_symbol(long long a, long long p) {
    if ((a % p) == 0) {
        return 0;
    }
    long long exponent = (p - 1) / 2;
    long long res = power(a, exponent, p);

    if (res == p - 1) {
        return -1;
    }
    // By Euler's Criterion, if a is not a multiple of p, the result
    // of the exponentiation must be 1 or p-1. So if not p-1, it is 1.
    return 1;
}


// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_pairs> <modulus> <seed>\n", argv[0]);
        exit(1);
    }

    num_pairs = atoi(argv[1]);
    modulus = atoll(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (num_pairs <= 0 || modulus <= 2 || (modulus % 2 == 0)) {
        fprintf(stderr, "Error: num_pairs must be > 0 and modulus must be an odd prime > 2.\n");
        exit(1);
    }

    mt_seed(seed);

    a_values = (long long *)malloc(num_pairs * sizeof(long long));
    if (a_values == NULL) {
        fprintf(stderr, "Failed to allocate memory for a_values.\n");
        exit(1);
    }

    l_symbols = (int *)malloc(num_pairs * sizeof(int));
    if (l_symbols == NULL) {
        fprintf(stderr, "Failed to allocate memory for l_symbols.\n");
        free(a_values);
        exit(1);
    }

    for (int i = 0; i < num_pairs; i++) {
        // Generate values 'a' in the range [0, modulus-1]
        a_values[i] = mt_rand() % modulus;
    }
    
    final_result = 0;
}

void run_computation() {
    long long sum = 0;
    
    // Compute all Legendre symbols
    for (int i = 0; i < num_pairs; i++) {
        l_symbols[i] = legendre_symbol(a_values[i], modulus);
    }

    // Accumulate results to prevent dead-code elimination
    for (int i = 0; i < num_pairs; i++) {
        sum += l_symbols[i];
    }
    final_result = sum;
}

void cleanup() {
    free(a_values);
    free(l_symbols);
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
    
    // Print the final accumulated result to stdout
    printf("%lld\n", final_result);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
