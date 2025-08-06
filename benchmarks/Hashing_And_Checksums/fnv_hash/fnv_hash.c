#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

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

// --- Benchmark Globals ---
unsigned char *data_buffer;
size_t data_buffer_size;
uint64_t final_result; // Use uint64_t to prevent overflow of summed hashes

// --- FNV-1a Hash Function (Helper) ---
// FNV-1a 32-bit constants
#define FNV_PRIME_32 16777619u
#define FNV_OFFSET_BASIS_32 2166136261u

// Calculates the FNV-1a hash of a given data block.
static inline uint32_t fnv1a_hash(const unsigned char *data, size_t len) {
    uint32_t hash = FNV_OFFSET_BASIS_32;
    for (size_t i = 0; i < len; ++i) {
        hash ^= data[i];
        hash *= FNV_PRIME_32;
    }
    return hash;
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <data_size_kb> <seed>\n", argv[0]);
        exit(1);
    }

    long long data_size_kb = atoll(argv[1]);
    uint32_t seed = (uint32_t)strtoul(argv[2], NULL, 10);

    if (data_size_kb <= 0) {
        fprintf(stderr, "FATAL: Invalid data_size_kb: %lld\n", data_size_kb);
        exit(1);
    }
    data_buffer_size = data_size_kb * 1024;

    data_buffer = (unsigned char *)malloc(data_buffer_size);
    if (data_buffer == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed for %zu bytes.\n", data_buffer_size);
        exit(1);
    }

    mt_seed(seed);

    // Populate the buffer with random bytes
    uint32_t r_val = 0;
    for (size_t i = 0; i < data_buffer_size; ++i) {
        if (i % 4 == 0) {
            r_val = mt_rand();
        }
        data_buffer[i] = (r_val >> ((i % 4) * 8)) & 0xFF;
    }
}

void run_computation() {
    const size_t CHUNK_SIZE = 256;
    size_t num_chunks = data_buffer_size / CHUNK_SIZE;
    size_t remainder = data_buffer_size % CHUNK_SIZE;
    uint64_t accumulator = 0;

    // Process full chunks
    for (size_t i = 0; i < num_chunks; ++i) {
        const unsigned char* chunk_start = data_buffer + (i * CHUNK_SIZE);
        accumulator += fnv1a_hash(chunk_start, CHUNK_SIZE);
    }

    // Process the remaining part if it exists
    if (remainder > 0) {
        const unsigned char* remainder_start = data_buffer + (num_chunks * CHUNK_SIZE);
        accumulator += fnv1a_hash(remainder_start, remainder);
    }
    
    final_result = accumulator;
}

void cleanup() {
    free(data_buffer);
    data_buffer = NULL;
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final accumulated hash to stdout
    printf("%llu\n", (unsigned long long)final_result);

    // Print the execution time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
