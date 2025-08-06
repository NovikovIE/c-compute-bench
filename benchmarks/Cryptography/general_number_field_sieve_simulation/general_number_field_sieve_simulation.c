#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

// MERSENNE TWISTER (DO NOT MODIFY)
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
// END MERSENNE TWISTER

// Benchmark-specific global data
int NUMBER_BIT_LENGTH;
unsigned int* factor_base = NULL;
unsigned char* sieve_array = NULL;
size_t factor_base_size;
size_t sieve_array_size;
long long final_result = 0;

// Simulates generating a factor base and sieving interval
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <number_bit_length> <seed>\n", argv[0]);
        exit(1);
    }

    NUMBER_BIT_LENGTH = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if (NUMBER_BIT_LENGTH <= 0) {
        fprintf(stderr, "FATAL: number_bit_length must be positive.\n");
        exit(1);
    }
    
    mt_seed(seed);

    // Determine sieve and factor base sizes based on the input parameter.
    // These constants are tuned to provide meaningful work for a modern CPU.
    sieve_array_size = (size_t)NUMBER_BIT_LENGTH * 2000000;
    factor_base_size = (size_t)NUMBER_BIT_LENGTH * 250;

    // Heuristic for Sieve of Eratosthenes limit to find enough primes.
    // A simple upper bound of N*15 is used for robustness.
    size_t prime_limit = factor_base_size * 15;
    if (factor_base_size > 0 && prime_limit / factor_base_size != 15) { // Overflow check
        prime_limit = SIZE_MAX;
    }
    if (prime_limit < 100) prime_limit = 100;

    // Generate primes using Sieve of Eratosthenes
    char* is_prime = (char*)malloc(prime_limit + 1);
    if (!is_prime) {
        fprintf(stderr, "FATAL: Failed to allocate memory for prime generation.\n");
        exit(1);
    }
    memset(is_prime, 1, prime_limit + 1);
    is_prime[0] = is_prime[1] = 0;
    for (size_t p = 2; p * p <= prime_limit; ++p) {
        if (is_prime[p]) {
            for (size_t i = p * p; i <= prime_limit; i += p) {
                is_prime[i] = 0;
            }
        }
    }

    // Allocate and populate the factor base
    factor_base = (unsigned int*)malloc(factor_base_size * sizeof(unsigned int));
    if (!factor_base) {
        fprintf(stderr, "FATAL: Failed to allocate memory for factor_base.\n");
        free(is_prime);
        exit(1);
    }

    size_t count = 0;
    for (size_t p = 2; p <= prime_limit && count < factor_base_size; ++p) {
        if (is_prime[p]) {
            factor_base[count++] = (unsigned int)p;
        }
    }
    
    if (count < factor_base_size) {
        fprintf(stderr, "FATAL: Prime limit heuristic too low. Found %zu, need %zu.\n", count, factor_base_size);
        free(is_prime);
        free(factor_base);
        exit(1);
    }
    free(is_prime);

    // Allocate the main sieve array
    sieve_array = (unsigned char*)calloc(sieve_array_size, sizeof(unsigned char));
    if (!sieve_array) {
        fprintf(stderr, "FATAL: Failed to allocate memory for sieve_array.\n");
        free(factor_base);
        exit(1);
    }
}

// Simulates the sieving stage of the General Number Field Sieve (GNFS)
void run_computation() {
    // For each prime in our factor base, iterate over the sieve array and
    // mark all its multiples. "Marking" is an increment here, which simulates
    // the logarithmic additions in a real sieve. This part is memory-intensive.
    for (size_t i = 0; i < factor_base_size; ++i) {
        unsigned int p = factor_base[i];
        for (size_t j = p; j < sieve_array_size; j += p) {
            sieve_array[j]++;
        }
    }

    // To produce a final result and avoid dead code elimination,
    // scan the array and compute a sum of all sieve values.
    // This simulates searching for "smooth" relations.
    long long temp_sum = 0;
    for (size_t i = 0; i < sieve_array_size; ++i) {
        temp_sum += sieve_array[i];
    }
    final_result = temp_sum;
}

void cleanup() {
    free(factor_base);
    free(sieve_array);
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

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
