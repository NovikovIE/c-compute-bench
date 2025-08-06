#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// A. Mersenne Twister (MT19937) PRNG (DO NOT MODIFY)
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

// B. Benchmark-specific data and functions
unsigned char *data_buffer;
size_t data_size;
uint32_t final_result;

// CRC-16-CCITT lookup table
static uint16_t crc16_table[256];

// Generates the CRC16 lookup table.
static void generate_crc16_table() {
    const uint16_t polynomial = 0x1021;
    for (uint16_t i = 0; i < 256; i++) {
        uint16_t crc = i << 8;
        for (int j = 0; j < 8; j++) {
            if (crc & 0x8000) {
                crc = (crc << 1) ^ polynomial;
            } else {
                crc = crc << 1;
            }
        }
        crc16_table[i] = crc;
    }
}

// C. Benchmark implementation functions
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <data_size_kb> <seed>\n", argv[0]);
        exit(1);
    }

    size_t data_size_kb = atol(argv[1]);
    uint32_t seed = atoi(argv[2]);

    data_size = data_size_kb * 1024;

    mt_seed(seed);
    generate_crc16_table();

    data_buffer = (unsigned char *)malloc(data_size);
    if (data_buffer == NULL) {
        fprintf(stderr, "Memory allocation failed for data_buffer\n");
        exit(1);
    }

    for (size_t i = 0; i < data_size; i += 4) {
        uint32_t r = mt_rand();
        // Distribute the 4 bytes of the random number into the buffer
        data_buffer[i] = (unsigned char)((r >> 24) & 0xFF);
        if (i + 1 < data_size) data_buffer[i + 1] = (unsigned char)((r >> 16) & 0xFF);
        if (i + 2 < data_size) data_buffer[i + 2] = (unsigned char)((r >> 8) & 0xFF);
        if (i + 3 < data_size) data_buffer[i + 3] = (unsigned char)(r & 0xFF);
    }
}

void run_computation() {
    const size_t chunk_size = 4096; // Process in 4KB chunks
    size_t num_chunks = data_size / chunk_size;
    uint32_t accumulator = 0;

    for (size_t i = 0; i < num_chunks; ++i) {
        uint16_t crc = 0xFFFF; // Initial value for CRC-16-CCITT-FALSE
        const unsigned char* chunk_start = data_buffer + (i * chunk_size);

        for (size_t j = 0; j < chunk_size; ++j) {
            uint8_t index = (uint8_t)(crc >> 8) ^ chunk_start[j];
            crc = (crc << 8) ^ crc16_table[index];
        }
        accumulator += crc;
    }

    final_result = accumulator;
}

void cleanup() {
    free(data_buffer);
}

// D. Main function, timing, and output
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final accumulated result to stdout
    printf("%u\n", final_result);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
