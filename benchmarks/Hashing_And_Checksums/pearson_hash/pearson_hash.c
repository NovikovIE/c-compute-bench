#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- Mersenne Twister (MT19937) Generator --- Do Not Modify ---
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
// --- End of MT19937 --- 

// --- Benchmark Globals ---
static size_t DATA_SIZE_BYTES;
static unsigned char* data;
static uint8_t T[256]; // Pearson hash lookup table
static uint8_t final_hash;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <data_size_kb> <seed>\n", argv[0]);
        exit(1);
    }

    long data_size_kb = atol(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    DATA_SIZE_BYTES = data_size_kb * 1024;

    mt_seed(seed);

    // 1. Generate the random permutation table T for Pearson hashing
    for (int i = 0; i < 256; i++) {
        T[i] = (uint8_t)i;
    }
    // Fisher-Yates shuffle
    for (int i = 255; i > 0; i--) {
        int j = mt_rand() % (i + 1);
        uint8_t temp = T[i];
        T[i] = T[j];
        T[j] = temp;
    }

    // 2. Allocate and fill the main data buffer
    data = (unsigned char*)malloc(DATA_SIZE_BYTES);
    if (data == NULL) {
        fprintf(stderr, "Failed to allocate memory for data buffer.\n");
        exit(1);
    }

    for (size_t i = 0; i < DATA_SIZE_BYTES; i += 4) {
        uint32_t r = mt_rand();
        if (i < DATA_SIZE_BYTES) data[i] = r & 0xFF;
        if (i + 1 < DATA_SIZE_BYTES) data[i+1] = (r >> 8) & 0xFF;
        if (i + 2 < DATA_SIZE_BYTES) data[i+2] = (r >> 16) & 0xFF;
        if (i + 3 < DATA_SIZE_BYTES) data[i+3] = (r >> 24) & 0xFF;
    }
}

void run_computation() {
    uint8_t h = 0;
    for (size_t i = 0; i < DATA_SIZE_BYTES; i++) {
        h = T[h ^ data[i]];
    }
    final_hash = h;
}

void cleanup() {
    free(data);
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

    // Print result to stdout
    printf("%u\n", (unsigned int)final_hash);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
