#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
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

// --- xxHash32 Implementation ---
// Based on the public domain implementation of xxHash by Yann Collet
#define XXH_PRIME32_1 2654435761U
#define XXH_PRIME32_2 2246822519U
#define XXH_PRIME32_3 3266489917U
#define XXH_PRIME32_4 668265263U
#define XXH_PRIME32_5 374761393U

static inline uint32_t XXH_rotl32(uint32_t x, int r) {
    return (x << r) | (x >> (32 - r));
}

static inline uint32_t XXH_read32(const void* memptr) {
    uint32_t val;
    memcpy(&val, memptr, sizeof(val));
    return val;
}

static uint32_t XXH32(const void* input, size_t len, uint32_t seed) {
    const uint8_t* p = (const uint8_t*)input;
    const uint8_t* const bEnd = p + len;
    uint32_t h32;

    if (len >= 16) {
        const uint8_t* const limit = bEnd - 16;
        uint32_t v1 = seed + XXH_PRIME32_1 + XXH_PRIME32_2;
        uint32_t v2 = seed + XXH_PRIME32_2;
        uint32_t v3 = seed + 0;
        uint32_t v4 = seed - XXH_PRIME32_1;

        do {
            v1 += XXH_read32(p) * XXH_PRIME32_2; v1 = XXH_rotl32(v1, 13); v1 *= XXH_PRIME32_1; p += 4;
            v2 += XXH_read32(p) * XXH_PRIME32_2; v2 = XXH_rotl32(v2, 13); v2 *= XXH_PRIME32_1; p += 4;
            v3 += XXH_read32(p) * XXH_PRIME32_2; v3 = XXH_rotl32(v3, 13); v3 *= XXH_PRIME32_1; p += 4;
            v4 += XXH_read32(p) * XXH_PRIME32_2; v4 = XXH_rotl32(v4, 13); v4 *= XXH_PRIME32_1; p += 4;
        } while (p <= limit);

        h32 = XXH_rotl32(v1, 1) + XXH_rotl32(v2, 7) + XXH_rotl32(v3, 12) + XXH_rotl32(v4, 18);
    } else {
        h32 = seed + XXH_PRIME32_5;
    }

    h32 += (uint32_t)len;

    while (p + 4 <= bEnd) {
        h32 += XXH_read32(p) * XXH_PRIME32_3; h32 = XXH_rotl32(h32, 17) * XXH_PRIME32_4; p += 4;
    }

    while (p < bEnd) {
        h32 += (*p) * XXH_PRIME32_5; h32 = XXH_rotl32(h32, 11) * XXH_PRIME32_1; p++;
    }

    h32 ^= h32 >> 15; h32 *= XXH_PRIME32_2; h32 ^= h32 >> 13; h32 *= XXH_PRIME32_3; h32 ^= h32 >> 16;

    return h32;
}

// --- Benchmark Globals ---
unsigned char *data_buffer;
size_t data_size_bytes;
uint32_t user_seed;
uint32_t final_hash;

// --- Benchmark Functions ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <data_size_kb> <seed>\n", argv[0]);
        exit(1);
    }

    unsigned long data_size_kb = strtoul(argv[1], NULL, 10);
    user_seed = (uint32_t)strtoul(argv[2], NULL, 10);

    data_size_bytes = data_size_kb * 1024;
    mt_seed(user_seed);

    data_buffer = (unsigned char *)malloc(data_size_bytes);
    if (!data_buffer) {
        fprintf(stderr, "FATAL: Memory allocation failed for data buffer of size %zu bytes.\n", data_size_bytes);
        exit(1);
    }
    
    unsigned char* p = data_buffer;
    for (size_t i = 0; i < data_size_bytes / 4; ++i) {
        uint32_t r = mt_rand();
        memcpy(p, &r, 4);
        p += 4;
    }
    size_t remainder = data_size_bytes % 4;
    if (remainder > 0) {
        uint32_t r = mt_rand();
        memcpy(p, &r, remainder);
    }
}

void run_computation() {
    final_hash = XXH32(data_buffer, data_size_bytes, user_seed);
}

void cleanup() {
    free(data_buffer);
    data_buffer = NULL;
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

    printf("%u\n", final_hash);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
