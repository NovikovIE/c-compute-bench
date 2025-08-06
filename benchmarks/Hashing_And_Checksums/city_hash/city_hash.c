#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

// --- Mersenne Twister (verbatim) ---
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
// --- End Mersenne Twister ---

// --- CityHash C Implementation --- 
// Relevant parts of CityHash (MIT License) for hashing 33-64 byte chunks.

static const uint64_t k0 = 0xc3a5c85c97cb3127ULL;
static const uint64_t k1 = 0xb492b66fbe98f273ULL;
static const uint64_t k2 = 0x9ae16a3b2f90404fULL;

static inline uint64_t Fetch64(const char *p) {
    uint64_t result;
    memcpy(&result, p, sizeof(result));
    return result;
}

static inline uint64_t Rotate(uint64_t val, int shift) {
    return shift == 0 ? val : ((val >> shift) | (val << (64 - shift)));
}

static inline uint64_t HashLen16(uint64_t u, uint64_t v, uint64_t mul) {
  uint64_t a = (u ^ v) * mul;
  a ^= (a >> 47);
  uint64_t b = (v ^ a) * mul;
  b ^= (b >> 47);
  b *= mul;
  return b;
}

static uint64_t CityHash64_33to64(const char *s, size_t len) {
  uint64_t mul = k2 + len * 2;
  uint64_t a = Fetch64(s) * k2;
  uint64_t b = Fetch64(s + 8);
  uint64_t c = Fetch64(s + len - 24);
  uint64_t d = Fetch64(s + len - 32);
  uint64_t e = Fetch64(s + 16) * k2;
  uint64_t f = Fetch64(s + 24) * 9;
  uint64_t g = Fetch64(s + len - 8);
  uint64_t h = Fetch64(s + len - 16) * mul;
  uint64_t u = Rotate(a + g, 43) + (Rotate(b, 30) + c) * 9;
  uint64_t v = (u + d)*mul + e;
  uint64_t w = Rotate(f+g, 42) + (u + v);
  uint64_t x = Rotate(e+h, 53) + w;
  uint64_t y = (x + d)*mul + f;
  uint64_t z = Rotate(h+a, 52) + y;
  return HashLen16(z, (z + e)*mul + g, mul);
}
// --- End CityHash C Implementation ---

// --- Benchmark Globals ---
#define NUM_HASHES 100000000
#define CHUNK_SIZE 64

typedef struct {
    char *data_buffer;
    size_t data_size;
    uint32_t *hash_indices;
    uint64_t final_hash;
} BenchmarkData;

static BenchmarkData g_data;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <data_size_kb> <seed>\n", argv[0]);
        exit(1);
    }

    int data_size_kb = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if (data_size_kb <= 0) {
        fprintf(stderr, "ERROR: data_size_kb must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    g_data.data_size = (size_t)data_size_kb * 1024;
    if (g_data.data_size <= CHUNK_SIZE) {
        fprintf(stderr, "ERROR: data_size must be greater than CHUNK_SIZE (%d bytes).\n", CHUNK_SIZE);
        exit(1);
    }
    
    g_data.data_buffer = (char*)malloc(g_data.data_size);
    if (!g_data.data_buffer) {
        perror("malloc data_buffer");
        exit(1);
    }

    for (size_t i = 0; i < g_data.data_size; i += sizeof(uint32_t)) {
        uint32_t r = mt_rand();
        size_t to_copy = (i + sizeof(uint32_t) <= g_data.data_size) ? sizeof(uint32_t) : g_data.data_size - i;
        memcpy(g_data.data_buffer + i, &r, to_copy);
    }

    g_data.hash_indices = (uint32_t*)malloc(NUM_HASHES * sizeof(uint32_t));
    if (!g_data.hash_indices) {
        perror("malloc hash_indices");
        free(g_data.data_buffer);
        exit(1);
    }
    
    uint32_t max_index = g_data.data_size - CHUNK_SIZE;
    for (int i = 0; i < NUM_HASHES; ++i) {
        g_data.hash_indices[i] = mt_rand() % max_index;
    }

    g_data.final_hash = 0;
}

void run_computation() {
    uint64_t accumulator = 0;
    for (int i = 0; i < NUM_HASHES; ++i) {
        uint32_t index = g_data.hash_indices[i];
        accumulator ^= CityHash64_33to64(g_data.data_buffer + index, CHUNK_SIZE);
    }
    g_data.final_hash = accumulator;
}

void cleanup() {
    free(g_data.data_buffer);
    free(g_data.hash_indices);
    g_data.data_buffer = NULL;
    g_data.hash_indices = NULL;
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%llu\n", g_data.final_hash);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}