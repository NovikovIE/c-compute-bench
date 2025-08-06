#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- MERSENNE TWISTER (DO NOT MODIFY) ---
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
// --- END MERSENNE TWISTER ---

// --- BENCHMARK-SPECIFIC STATE ---
static uint8_t *data_buffer = NULL;
static size_t data_size;
static uint8_t key[16];
static uint64_t final_result = 0;

// --- SIPHASH IMPLEMENTATION HELPERS ---
#define ROTL64(x, b) (uint64_t)(((x) << (b)) | ((x) >> (64 - (b))))

// Helper to read a little-endian uint64_t from a byte array
static inline uint64_t U8TO64_LE(const uint8_t *p) {
    return ((uint64_t)p[0]) | ((uint64_t)p[1] << 8) |
           ((uint64_t)p[2] << 16) | ((uint64_t)p[3] << 24) |
           ((uint64_t)p[4] << 32) | ((uint64_t)p[5] << 40) |
           ((uint64_t)p[6] << 48) | ((uint64_t)p[7] << 56);
}

#define SIPROUND                                     \
    do {                                             \
        v0 += v1; v1 = ROTL64(v1, 13); v1 ^= v0;      \
        v0 = ROTL64(v0, 32);                         \
        v2 += v3; v3 = ROTL64(v3, 16); v3 ^= v2;      \
        v0 += v3; v3 = ROTL64(v3, 21); v3 ^= v0;      \
        v2 += v1; v1 = ROTL64(v1, 17); v1 ^= v2;      \
        v2 = ROTL64(v2, 32);                         \
    } while (0)

// SipHash-2-4 implementation
static uint64_t do_siphash(const uint8_t *in, const size_t inlen, const uint8_t *k) {
    uint64_t v0 = 0x736f6d6570736575ULL;
    uint64_t v1 = 0x646f72616e646f6dULL;
    uint64_t v2 = 0x6c7967656e657261ULL;
    uint64_t v3 = 0x7465646279746573ULL;
    uint64_t k0 = U8TO64_LE(k);
    uint64_t k1 = U8TO64_LE(k + 8);
    uint64_t m;
    const uint8_t *end = in + inlen - (inlen % 8);
    const int left = inlen & 7;
    uint64_t b = ((uint64_t)inlen) << 56;

    v3 ^= k1; v2 ^= k0; v1 ^= k1; v0 ^= k0;

    for (; in != end; in += 8) {
        m = U8TO64_LE(in);
        v3 ^= m;
        SIPROUND; SIPROUND;
        v0 ^= m;
    }

    switch (left) {
        case 7: b |= ((uint64_t)in[6]) << 48;
        case 6: b |= ((uint64_t)in[5]) << 40;
        case 5: b |= ((uint64_t)in[4]) << 32;
        case 4: b |= ((uint64_t)in[3]) << 24;
        case 3: b |= ((uint64_t)in[2]) << 16;
        case 2: b |= ((uint64_t)in[1]) << 8;
        case 1: b |= ((uint64_t)in[0]); break;
        case 0: break;
    }

    v3 ^= b;
    SIPROUND; SIPROUND;
    v0 ^= b;

    v2 ^= 0xff;
    SIPROUND; SIPROUND; SIPROUND; SIPROUND;

    return v0 ^ v1 ^ v2 ^ v3;
}

// --- BENCHMARK ROUTINES ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <data_size_kb> <seed>\n", argv[0]);
        exit(1);
    }

    long data_size_kb = atol(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    data_size = (size_t)data_size_kb * 1024;
    mt_seed(seed);

    data_buffer = (uint8_t *)malloc(data_size);
    if (data_buffer == NULL) {
        fprintf(stderr, "Failed to allocate memory for data buffer.\n");
        exit(1);
    }

    // Populate data buffer with random bytes
    for (size_t i = 0; i < data_size / 4; ++i) {
        uint32_t random_val = mt_rand();
        ((uint32_t*)data_buffer)[i] = random_val;
    }
    
    // Generate 128-bit (16-byte) key from the PRNG
    for (int i = 0; i < 4; ++i) {
        uint32_t random_val = mt_rand();
        ((uint32_t*)key)[i] = random_val;
    }

    final_result = 0;
}

void run_computation() {
    const size_t CHUNK_SIZE = 256;
    uint64_t accumulated_hash = 0;
    size_t num_chunks = data_size / CHUNK_SIZE;

    for (size_t i = 0; i < num_chunks; ++i) {
        accumulated_hash ^= do_siphash(data_buffer + (i * CHUNK_SIZE), CHUNK_SIZE, key);
    }
    final_result = accumulated_hash;
}

void cleanup() {
    if (data_buffer) {
        free(data_buffer);
        data_buffer = NULL;
    }
}

// --- MAIN ---

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%llu\n", (unsigned long long)final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
