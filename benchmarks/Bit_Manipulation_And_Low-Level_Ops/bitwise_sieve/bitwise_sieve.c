#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

// --- Mersenne Twister (MT19937) Generator ---
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
// --- End of Mersenne Twister ---

// --- Benchmark Globals ---
typedef struct {
    uint32_t *sieve_array;  // Bit array for the sieve
    int upper_limit;        // The upper bound to search for primes
    long long prime_count;  // The final result
} BenchmarkData;

static BenchmarkData g_data;

// --- Benchmark Functions ---

/**
 * @brief Sets up the Sieve of Eratosthenes benchmark.
 *
 * Parses command-line arguments for the upper limit and a random seed.
 * Allocates a bit array on the heap large enough to hold primality
 * information for all numbers up to 'upper_limit'.
 * Although the sieve itself is deterministic, the random seed is initialized
 * as per the requirements.
 * The bit array is initialized to 0 using calloc, representing all numbers as potentially prime.
 * 0 and 1 are explicitly marked as not prime.
 */
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <upper_limit> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.upper_limit = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);
    
    if (g_data.upper_limit < 2) {
        fprintf(stderr, "Error: upper_limit must be at least 2.\n");
        exit(1);
    }

    // Seed the random number generator (as required, though not used for data generation)
    mt_seed(seed);

    // Calculate the size of the bit array. Each uint32_t can hold 32 flags.
    size_t num_words = (size_t)g_data.upper_limit / 32 + 1;
    g_data.sieve_array = (uint32_t*) calloc(num_words, sizeof(uint32_t));
    if (g_data.sieve_array == NULL) {
        fprintf(stderr, "Error: Memory allocation failed for sieve_array.\n");
        exit(1);
    }
    
    g_data.prime_count = 0;

    // Mark 0 and 1 as not prime.
    // Set bit 0:
    g_data.sieve_array[0] |= (1U << 0);
    // Set bit 1:
    g_data.sieve_array[0] |= (1U << 1);
}

/**
 * @brief Runs the core computation: the Sieve of Eratosthenes.
 *
 * This function identifies all prime numbers up to `upper_limit` using a
 * bitwise sieve. 
 * 1. It iterates from p = 2 up to sqrt(upper_limit).
 * 2. If a number p is found to be prime (its corresponding bit is 0),
 *    all multiples of p are marked as composite (their bits are set to 1).
 *    The marking starts from p*p, as smaller multiples would have been
 *    marked by smaller primes.
 * 3. After the sieve process, it counts all the numbers that remain marked
 *    as prime (bit is 0) and stores the total in `g_data.prime_count`.
 */
void run_computation() {
    int limit = g_data.upper_limit;
    uint32_t *sieve = g_data.sieve_array;
    int sqrt_limit = (int)sqrt((double)limit);

    // Mark multiples of primes
    for (int p = 2; p <= sqrt_limit; ++p) {
        // Check if p is prime (if its bit is 0)
        if (!((sieve[(size_t)p / 32] >> (p % 32)) & 1U)) {
            // Mark all multiples of p starting from p*p
            for (long long i = (long long)p * p; i <= limit; i += p) {
                sieve[(size_t)i / 32] |= (1U << (i % 32));
            }
        }
    }

    // Count the primes
    long long count = 0;
    for (int i = 2; i <= limit; ++i) {
        if (!((sieve[(size_t)i / 32] >> (i % 32)) & 1U)) {
            count++;
        }
    }
    g_data.prime_count = count;
}


/**
 * @brief Frees all memory allocated during setup.
 */
void cleanup() {
    free(g_data.sieve_array);
    g_data.sieve_array = NULL;
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

    // Print the final result (total number of primes) to stdout
    printf("%lld\n", g_data.prime_count);

    // Print the execution time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
