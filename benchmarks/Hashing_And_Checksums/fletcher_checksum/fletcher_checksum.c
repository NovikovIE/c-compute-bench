#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- BEGIN MERSENNE TWISTER (MT19937) ---
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

// Global variables for benchmark data
static uint8_t *data;
static size_t data_size_bytes;
static uint16_t final_checksum;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <data_size_kb> <seed>\n", argv[0]);
        exit(1);
    }

    long data_size_kb = atol(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    data_size_bytes = data_size_kb * 1024;
    if (data_size_bytes <= 0) {
        fprintf(stderr, "FATAL: Invalid data size.\n");
        exit(1);
    }

    mt_seed(seed);

    data = (uint8_t *)malloc(data_size_bytes);
    if (data == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    for (size_t i = 0; i < data_size_bytes; ++i) {
        data[i] = mt_rand() & 0xFF;
    }
}

void run_computation() {
    // Implementation of Fletcher-16 checksum with 32-bit accumulators
    // to reduce the frequency of modulo operations for performance.
    uint32_t c0 = 0;
    uint32_t c1 = 0;
    
    // Using a large block size is possible with 32-bit accumulators.
    // The largest k such that k additions of 255 do not overflow sum1 (c0) is
    // k*255 < 2^32, k < 16,843,009.
    // The largest k for sum2 (c1) is roughly k^2*127.5 < 2^32, so k < 5804.
    // We use a safe block size of 5800.
    const size_t block_size = 5800;
    size_t remaining_bytes = data_size_bytes;
    size_t i = 0;

    while (remaining_bytes > 0) {
        size_t current_block_len = remaining_bytes > block_size ? block_size : remaining_bytes;
        remaining_bytes -= current_block_len;
        
        for (size_t j = 0; j < current_block_len; ++j) {
            c0 += data[i++];
            c1 += c0;
        }
        
        c0 %= 255;
        c1 %= 255;
    }

    final_checksum = (uint16_t)((c1 << 8) | c0);
}

void cleanup() {
    if (data) {
        free(data);
        data = NULL;
    }
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final checksum to stdout
    printf("%u\n", (unsigned int)final_checksum);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
