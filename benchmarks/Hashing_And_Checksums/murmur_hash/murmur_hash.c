#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- BEGIN: Mersenne Twister (Do Not Modify) ---
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
// --- END: Mersenne Twister ---

// --- Benchmark Globals ---
unsigned char *data;
size_t data_size_bytes;
uint32_t final_result;

// --- MurmurHash3 Helper ---
// 32-bit Murmur3 hash function implementation
static uint32_t murmur3_32(const uint8_t* key, size_t len, uint32_t seed) {
    uint32_t h = seed;
    const uint32_t c1 = 0xcc9e2d51;
    const uint32_t c2 = 0x1b873593;
    const int r1 = 15;
    const int r2 = 13;
    const uint32_t m = 5;
    const uint32_t n = 0xe6546b64;

    const int nblocks = len / 4;
    const uint32_t *blocks = (const uint32_t *)key;
    int i;
    for (i = 0; i < nblocks; i++) {
        uint32_t k = blocks[i];
        k *= c1;
        k = (k << r1) | (k >> (32 - r1));
        k *= c2;
        h ^= k;
        h = (h << r2) | (h >> (32 - r2));
        h = h * m + n;
    }

    const uint8_t *tail = (const uint8_t *)(key + nblocks * 4);
    uint32_t k1 = 0;

    switch (len & 3) {
        case 3: k1 ^= tail[2] << 16;
        case 2: k1 ^= tail[1] << 8;
        case 1: k1 ^= tail[0];
                k1 *= c1;
                k1 = (k1 << r1) | (k1 >> (32 - r1));
                k1 *= c2;
                h ^= k1;
    }

    h ^= len;
    h ^= h >> 16;
    h *= 0x85ebca6b;
    h ^= h >> 13;
    h *= 0xc2b2ae35;
    h ^= h >> 16;

    return h;
}


void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <data_size_kb> <seed>\n", argv[0]);
        exit(1);
    }

    long data_size_kb = atol(argv[1]);
    uint32_t seed = atoi(argv[2]);

    data_size_bytes = (size_t)data_size_kb * 1024;
    data = (unsigned char *)malloc(data_size_bytes);
    if (data == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed for data buffer.\n");
        exit(1);
    }

    mt_seed(seed);

    // Populate the data buffer with random bytes
    for (size_t i = 0; i < data_size_bytes; i++) {
        data[i] = (unsigned char)(mt_rand() & 0xFF);
    }
}

void run_computation() {
    const size_t CHUNK_SIZE = 65536; // 64KB chunks
    const uint32_t HASH_SEED = 0x9747b28c; // A fixed seed for the hash function
    uint32_t accumulated_hash = 0;
    size_t offset = 0;

    while (offset < data_size_bytes) {
        size_t current_chunk_size = (data_size_bytes - offset < CHUNK_SIZE) ? (data_size_bytes - offset) : CHUNK_SIZE;
        uint32_t chunk_hash = murmur3_32(data + offset, current_chunk_size, HASH_SEED);
        accumulated_hash ^= chunk_hash; // Accumulate results to prevent dead code elimination
        offset += current_chunk_size;
    }

    final_result = accumulated_hash;
}

void cleanup() {
    free(data);
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
    printf("%u\n", final_result);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
