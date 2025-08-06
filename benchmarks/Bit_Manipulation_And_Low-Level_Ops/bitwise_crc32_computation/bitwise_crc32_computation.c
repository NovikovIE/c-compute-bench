#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

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
// --- END Mersenne Twister ---

// --- Benchmark Globals ---
size_t data_size_bytes;
unsigned char* data_buffer;
uint32_t final_crc_result;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <data_size_bytes> <seed>\n", argv[0]);
        exit(1);
    }

    data_size_bytes = atol(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if (data_size_bytes <= 0) {
        fprintf(stderr, "FATAL: data_size_bytes must be a positive integer.\n");
        exit(1);
    }
    
    mt_seed(seed);

    data_buffer = (unsigned char*)malloc(data_size_bytes);
    if (data_buffer == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed for data_buffer.\n");
        exit(1);
    }

    for (size_t i = 0; i < data_size_bytes; i++) {
        data_buffer[i] = (unsigned char)(mt_rand() & 0xFF);
    }
}

void run_computation() {
    // This is a naive, bit-by-bit implementation of CRC32.
    // A slow implementation is chosen to emphasize bitwise operations.
    // Polynomial: 0x04C11DB7, Reflected: 0xEDB88320
    const uint32_t polynomial = 0xEDB88320;
    uint32_t crc = 0xFFFFFFFF; // Initial value

    for (size_t i = 0; i < data_size_bytes; i++) {
        crc ^= data_buffer[i];
        for (int j = 0; j < 8; j++) {
            if (crc & 1) {
                crc = (crc >> 1) ^ polynomial;
            } else {
                crc >>= 1;
            }
        }
    }

    final_crc_result = ~crc; // Final XOR
}

void cleanup() {
    free(data_buffer);
}

// --- Main ---

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
    // Print result to stdout
    printf("%u\n", final_crc_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
