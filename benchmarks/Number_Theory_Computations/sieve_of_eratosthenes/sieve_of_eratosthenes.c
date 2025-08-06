#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

// --- Mersenne Twister (DO NOT MODIFY) ---
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
// --- End of Mersenne Twister ---

// --- Benchmark Globals ---
static unsigned int N;
static bool *is_prime;
static long long prime_sum;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <upper_limit> <seed>\n", argv[0]);
        exit(1);
    }

    N = (unsigned int)atol(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);
    mt_seed(seed);

    if (N < 2) {
        fprintf(stderr, "Error: upper_limit must be at least 2.\n");
        exit(1);
    }

    is_prime = (bool *)malloc((N + 1) * sizeof(bool));
    if (is_prime == NULL) {
        fprintf(stderr, "Error: Failed to allocate memory for the sieve.\n");
        exit(1);
    }

    // The Sieve of Eratosthenes algorithm requires an initial assumption that all numbers are prime.
    // This is part of the data setup.
    for (unsigned int i = 0; i <= N; i++) {
        is_prime[i] = true;
    }
}

void run_computation() {
    is_prime[0] = is_prime[1] = false;
    
    for (unsigned long long p = 2; p * p <= N; p++) {
        if (is_prime[p]) {
            for (unsigned long long i = p * p; i <= N; i += p) {
                is_prime[i] = false;
            }
        }
    }

    // To prevent dead code elimination, we sum the primes found.
    prime_sum = 0;
    for (unsigned int i = 2; i <= N; i++) {
        if (is_prime[i]) {
            prime_sum += i;
        }
    }
}

void cleanup() {
    free(is_prime);
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

    // Print result to stdout
    printf("%lld\n", prime_sum);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
