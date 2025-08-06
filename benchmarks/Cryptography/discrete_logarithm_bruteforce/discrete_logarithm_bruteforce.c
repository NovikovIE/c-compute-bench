#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- Mersenne Twister (MT19937) Generator (Verbatim) ---
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
// Problem parameters: find x such that g^x ≡ h (mod p)
unsigned long long p; // Prime modulus
unsigned long long g; // Generator (base)
unsigned long long h; // Target value

// Result of the computation
unsigned long long found_exponent;

// --- Utility Functions ---

// Computes (base^exp) % modulus using modular exponentiation.
// Uses __int128 to prevent overflow during intermediate multiplication.
static unsigned long long power(unsigned long long base, unsigned long long exp, unsigned long long modulus) {
    unsigned long long res = 1;
    base %= modulus;
    while (exp > 0) {
        if (exp % 2 == 1) res = ((__int128)res * base) % modulus;
        base = ((__int128)base * base) % modulus;
        exp /= 2;
    }
    return res;
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <prime_field_bit_length> <seed>\n", argv[0]);
        exit(1);
    }

    int exponent_bit_length = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);
    mt_seed(seed);

    if (exponent_bit_length < 1 || exponent_bit_length > 30) {
        fprintf(stderr, "FATAL: prime_field_bit_length (interpreted as exponent bits) must be between 1 and 30.\n");
        exit(1);
    }

    // For this benchmark, we fix the prime field and generator to keep the
    // modular arithmetic cost constant. The 'prime_field_bit_length' parameter
    // controls the size of the secret exponent, and thus the brute-force search space.
    
    // p is the Mersenne prime 2^31 - 1, which fits in a 64-bit unsigned integer.
    p = 2147483647;
    g = 3; // A common base/generator

    // Determine the secret exponent 'x' we need to find.
    // Its size is determined by the input parameter.
    unsigned long long max_exponent = 1ULL << exponent_bit_length;
    // We want a non-trivial exponent (>= 2) to search for.
    unsigned long long secret_x = (mt_rand() % (max_exponent - 2)) + 2;

    // Calculate h = g^x (mod p). This is the value we give to the computation.
    h = power(g, secret_x, p);

    // Initialize the result variable.
    found_exponent = 0; // 0 indicates not found
}

void run_computation() {
    // Brute-force search for the exponent x by trial multiplication.
    // We start with x=1 and compute g^x iteratively until we find a match with h.
    unsigned long long current_power = g; // Starts at g^1
    for (unsigned long long x = 1; x < p; ++x) {
        if (current_power == h) {
            found_exponent = x;
            return;
        }
        // Calculate next power: g^(x+1) = (g^x * g) mod p
        current_power = ((__int128)current_power * g) % p;
    }
}

void cleanup() {
    // No heap-allocated memory to free in this benchmark.
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

    // Print the found exponent to stdout
    printf("%llu\n", found_exponent);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
