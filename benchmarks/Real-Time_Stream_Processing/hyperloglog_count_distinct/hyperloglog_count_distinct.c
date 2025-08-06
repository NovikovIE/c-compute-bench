#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <string.h>

// --- Mersenne Twister (MT19937) Generator (Provided) ---
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
// --- End of MT19937 ---

// --- HyperLogLog Implementation Details ---
// HLL uses a probabilistic approach to estimate the number of distinct elements.
// P defines the precision. The number of registers 'm' will be 2^P.
// Higher P means more memory and better accuracy.
// P=14 means m=16384 registers, a common setting.
#define HLL_P 14
#define HLL_M (1 << HLL_P)

// Correction constant alpha, depends on m.
#define HLL_ALPHA (0.7213 / (1.0 + 1.079 / HLL_M))

// Use built-in for leading zero count on GCC/Clang for performance
#if defined(__GNUC__) || defined(__clang__)
#define COUNT_LEADING_ZEROS(x) ((x) == 0 ? 32 : __builtin_clz(x))
#else
// Fallback for other compilers
static int count_leading_zeros_fallback(uint32_t x) {
    if (x == 0) return 32;
    int n = 1;
    if ((x >> 16) == 0) { n += 16; x <<= 16; }
    if ((x >> 24) == 0) { n += 8;  x <<= 8;  }
    if ((x >> 28) == 0) { n += 4;  x <<= 4;  }
    if ((x >> 30) == 0) { n += 2;  x <<= 2;  }
    n -= (x >> 31);
    return n;
}
#define COUNT_LEADING_ZEROS(x) count_leading_zeros_fallback(x)
#endif

// Simple integer hash function (Murmur3 finalizer)
static inline uint32_t hash_u32(uint32_t k) {
    k ^= k >> 16;
    k *= 0x85ebca6b;
    k ^= k >> 13;
    k *= 0xc2b2ae35;
    k ^= k >> 16;
    return k;
}

// --- Benchmark Globals and Functions ---

typedef struct {
    long total_events;
    long num_unique_keys;
    unsigned int seed;

    // Input data stream simulation
    uint32_t *events;

    // HyperLogLog state
    uint8_t *registers;

    // Final result
    long estimated_cardinality;
} BenchmarkData;

static BenchmarkData g_data;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s total_events num_unique_keys seed\n", argv[0]);
        exit(1);
    }

    g_data.total_events = atol(argv[1]);
    g_data.num_unique_keys = atol(argv[2]);
    g_data.seed = atoi(argv[3]);

    if (g_data.total_events <= 0 || g_data.num_unique_keys <= 0) {
        fprintf(stderr, "FATAL: total_events and num_unique_keys must be positive.\n");
        exit(1);
    }

    mt_seed(g_data.seed);

    // Allocate memory for the event stream
    g_data.events = (uint32_t *)malloc(g_data.total_events * sizeof(uint32_t));
    if (g_data.events == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed for events.\n");
        exit(1);
    }

    // Generate the event stream with keys up to num_unique_keys
    for (long i = 0; i < g_data.total_events; ++i) {
        g_data.events[i] = mt_rand() % g_data.num_unique_keys;
    }

    // Allocate and initialize HLL registers
    g_data.registers = (uint8_t *)malloc(HLL_M * sizeof(uint8_t));
    if (g_data.registers == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed for registers.\n");
        free(g_data.events);
        exit(1);
    }
    memset(g_data.registers, 0, HLL_M * sizeof(uint8_t));
}

void run_computation() {
    // --- Part 1: Process the stream and update HLL registers ---
    for (long i = 0; i < g_data.total_events; ++i) {
        uint32_t item = g_data.events[i];
        uint32_t h = hash_u32(item);

        // Use first P bits for register index
        uint32_t index = h >> (32 - HLL_P);
        // The rest of the hash for counting leading zeros
        uint32_t remainder = h << HLL_P;

        // Count leading zeros (ρ(w)) + 1
        uint8_t p = COUNT_LEADING_ZEROS(remainder) + 1;

        // Update register if new max is found
        if (p > g_data.registers[index]) {
            g_data.registers[index] = p;
        }
    }

    // --- Part 2: Calculate final cardinality estimate ---
    double sum_inv_pow2 = 0.0;
    int zero_registers = 0;
    for (int j = 0; j < HLL_M; j++) {
        if (g_data.registers[j] == 0) {
            zero_registers++;
            sum_inv_pow2 += 1.0;
        } else {
            sum_inv_pow2 += 1.0 / (1ULL << g_data.registers[j]);
        }
    }

    double estimate = HLL_ALPHA * HLL_M * HLL_M / sum_inv_pow2;

    // Apply low-cardinality correction if needed (LinearCounting)
    if (estimate <= 2.5 * HLL_M && zero_registers > 0) {
        estimate = HLL_M * log((double)HLL_M / zero_registers);
    }
    // Note: A full implementation would also have a high-cardinality correction.
    // For this benchmark, we omit it as the core computation is what we are measuring.

    g_data.estimated_cardinality = (long)(estimate + 0.5);
}

void cleanup() {
    free(g_data.events);
    free(g_data.registers);
    g_data.events = NULL;
    g_data.registers = NULL;
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%ld\n", g_data.estimated_cardinality);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
