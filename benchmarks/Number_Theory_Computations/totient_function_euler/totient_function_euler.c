#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

// Mersenne Twister (MT19937) Generator
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

// Benchmark-specific global data
static int num_samples;
static unsigned int upper_limit;
static unsigned int* numbers_to_test;
static long long final_result;

// Function to calculate Euler's totient function φ(n)
// φ(n) = n * Π(1 - 1/p) for all distinct prime factors p of n
long long euler_totient(unsigned int n) {
    if (n == 0) return 0;
    long long result = n;
    unsigned int p;

    // Handle factor 2 separately for efficiency
    if (n % 2 == 0) {
        result -= result / 2;
        while (n % 2 == 0) {
            n /= 2;
        }
    }

    // Check for odd factors from 3 up to sqrt(n)
    for (p = 3; p * p <= n; p += 2) {
        if (n % p == 0) {
            result -= result / p;
            while (n % p == 0) {
                n /= p;
            }
        }
    }

    // If n is still greater than 1, it must be a prime factor itself
    if (n > 1) {
        result -= result / n;
    }
    return result;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_samples> <upper_limit> <seed>\n", argv[0]);
        exit(1);
    }

    num_samples = atoi(argv[1]);
    upper_limit = (unsigned int)atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (num_samples <= 0 || upper_limit < 2) {
        fprintf(stderr, "Error: num_samples > 0 and upper_limit >= 2 required.\n");
        exit(1);
    }

    mt_seed(seed);

    numbers_to_test = (unsigned int*)malloc(num_samples * sizeof(unsigned int));
    if (numbers_to_test == NULL) {
        fprintf(stderr, "Error: Failed to allocate memory for test numbers.\n");
        exit(1);
    }

    for (int i = 0; i < num_samples; i++) {
        // Generate numbers in the range [2, upper_limit]
        numbers_to_test[i] = (mt_rand() % (upper_limit - 1)) + 2;
    }

    final_result = 0;
}

void run_computation() {
    long long sum = 0;
    for (int i = 0; i < num_samples; i++) {
        sum += euler_totient(numbers_to_test[i]);
    }
    final_result = sum;
}

void cleanup() {
    free(numbers_to_test);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout
    printf("%lld\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}