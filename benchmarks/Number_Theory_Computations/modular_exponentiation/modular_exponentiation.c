#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- Mersenne Twister (MT19937) ---
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
    uint64_t base;
    uint64_t exponent;
    uint64_t modulus;
} ModExpInput;

static ModExpInput* g_inputs = NULL;
static int g_num_calculations;
static uint64_t g_final_result = 0;

// --- Helper Functions ---

// Generate a random 64-bit integer
uint64_t mt_rand64() {
    uint64_t high = mt_rand();
    uint64_t low = mt_rand();
    return (high << 32) | low;
}

// Modular exponentiation (base^exp) % mod using binary exponentiation.
// Uses 128-bit integers for intermediate products to prevent overflow,
// which is required for 64-bit inputs.
uint64_t modular_pow(uint64_t base, uint64_t exp, uint64_t mod) {
    unsigned __int128 result = 1;
    unsigned __int128 b = base;

    if (mod == 1) return 0; // The only result mod 1 is 0

    b %= mod;
    while (exp > 0) {
        if (exp % 2 == 1) {
            result = (result * b) % mod;
        }
        b = (b * b) % mod;
        exp /= 2;
    }
    return (uint64_t)result;
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_calculations> <seed>\n", argv[0]);
        exit(1);
    }

    g_num_calculations = atoi(argv[1]);
    uint32_t seed = (uint32_t)strtoul(argv[2], NULL, 10);

    if (g_num_calculations <= 0) {
        fprintf(stderr, "ERROR: Number of calculations must be a positive integer.\n");
        exit(1);
    }

    mt_seed(seed);

    g_inputs = (ModExpInput*)malloc(g_num_calculations * sizeof(ModExpInput));
    if (g_inputs == NULL) {
        fprintf(stderr, "FATAL: Failed to allocate memory for inputs.\n");
        exit(1);
    }

    for (int i = 0; i < g_num_calculations; ++i) {
        g_inputs[i].base = mt_rand64();
        g_inputs[i].exponent = mt_rand64();
        
        // Modulus must be > 1 for the operation to be well-defined.
        do {
            g_inputs[i].modulus = mt_rand64();
        } while (g_inputs[i].modulus <= 1);
    }
}

void run_computation() {
    uint64_t accumulator = 0;
    for (int i = 0; i < g_num_calculations; ++i) {
        uint64_t result = modular_pow(
            g_inputs[i].base,
            g_inputs[i].exponent,
            g_inputs[i].modulus
        );
        // XOR result into accumulator to prevent dead-code elimination.
        accumulator ^= result;
    }
    g_final_result = accumulator;
}

void cleanup() {
    free(g_inputs);
    g_inputs = NULL;
}

// --- Main --- 

int main(int argc, char* argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%llu\n", (unsigned long long)g_final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
