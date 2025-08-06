#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <string.h>

// START: Mersenne Twister (Do Not Modify)
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
// END: Mersenne Twister

// --- Benchmark Globals ---
unsigned char* g_compressed_data = NULL;
unsigned char* g_decompressed_data = NULL;
size_t g_compressed_size = 0;
size_t g_decompressed_size = 0;
unsigned int g_final_checksum = 0;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <input_size_mb> <run_pattern_factor> <seed>\n", argv[0]);
        exit(1);
    }

    long input_size_mb = atol(argv[1]);
    int run_pattern_factor = atoi(argv[2]);
    uint32_t seed = (uint32_t)strtoul(argv[3], NULL, 10);

    if (input_size_mb <= 0 || run_pattern_factor <= 0) {
        fprintf(stderr, "FATAL: input_size_mb and run_pattern_factor must be positive.\n");
        exit(1);
    }
    if (run_pattern_factor > 255) {
        fprintf(stderr, "WARN: run_pattern_factor > 255 is capped at 255 for run length generation.\n");
        run_pattern_factor = 255;
    }

    mt_seed(seed);

    g_decompressed_size = (size_t)input_size_mb * 1024 * 1024;

    // Worst case for RLE is 2 bytes (count, value) for every 1 byte of original data.
    // Allocate a buffer guaranteed to be large enough.
    size_t max_compressed_size = g_decompressed_size * 2;

    g_decompressed_data = (unsigned char*)malloc(g_decompressed_size);
    g_compressed_data = (unsigned char*)malloc(max_compressed_size);

    if (!g_decompressed_data || !g_compressed_data) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Generate the RLE-compressed data stream
    size_t decompressed_bytes_generated = 0;
    size_t compressed_idx = 0;

    while (decompressed_bytes_generated < g_decompressed_size) {
        // Generate a run length from 1 to run_pattern_factor.
        unsigned int run_length = 1 + (mt_rand() % run_pattern_factor);

        // Ensure this run does not overshoot the total decompressed size.
        if (decompressed_bytes_generated + run_length > g_decompressed_size) {
            run_length = g_decompressed_size - decompressed_bytes_generated;
        }

        if (run_length == 0) {
            break; // Finished generating
        }

        unsigned char value = (unsigned char)(mt_rand() & 0xFF);
        
        g_compressed_data[compressed_idx++] = (unsigned char)run_length;
        g_compressed_data[compressed_idx++] = value;

        decompressed_bytes_generated += run_length;
    }

    g_compressed_size = compressed_idx;
}

void run_computation() {
    size_t src_idx = 0;
    size_t dest_idx = 0;

    while (src_idx < g_compressed_size) {
        unsigned char run_length = g_compressed_data[src_idx++];
        unsigned char value = g_compressed_data[src_idx++];
        
        // Decompress the run. memset is highly optimized for this task.
        // We assume setup generated valid data, so no overflow check is needed here for speed.
        memset(g_decompressed_data + dest_idx, value, run_length);
        dest_idx += run_length;
    }
    
    // To prevent dead code elimination by the compiler, we calculate a checksum
    // of the entire decompressed data block.
    unsigned int checksum = 0;
    for (size_t i = 0; i < g_decompressed_size; ++i) {
        checksum += g_decompressed_data[i];
    }
    g_final_checksum = checksum;
}

void cleanup() {
    free(g_compressed_data);
    g_compressed_data = NULL;
    free(g_decompressed_data);
    g_decompressed_data = NULL;
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final result to stdout
    printf("%u\n", g_final_checksum);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
