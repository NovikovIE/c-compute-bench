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
// --- END MERSENNE TWISTER ---

// --- GLOBAL BENCHMARK DATA ---
static uint8_t *g_input_data = NULL;
static uint8_t *g_output_data = NULL;
static size_t g_num_elements = 0;
static size_t g_total_bytes = 0;
static int g_word_size_bytes = 0;
static long long g_final_checksum = 0;
// --- END GLOBAL DATA ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <input_size_mb> <word_size_bytes> <seed>\n", argv[0]);
        exit(1);
    }

    size_t input_size_mb = atol(argv[1]);
    g_word_size_bytes = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    if (g_word_size_bytes != 1 && g_word_size_bytes != 2 && g_word_size_bytes != 4 && g_word_size_bytes != 8) {
        fprintf(stderr, "FATAL: word_size_bytes must be 1, 2, 4, or 8.\n");
        exit(1);
    }

    g_total_bytes = input_size_mb * 1024 * 1024;
    // Align total_bytes to be a multiple of word_size_bytes
    g_total_bytes -= g_total_bytes % g_word_size_bytes;
    g_num_elements = g_total_bytes / g_word_size_bytes;

    if (g_num_elements == 0) {
        fprintf(stderr, "FATAL: Calculated 0 elements. Input size might be too small for the word size.\n");
        exit(1);
    }

    g_input_data = (uint8_t*)malloc(g_total_bytes);
    g_output_data = (uint8_t*)malloc(g_total_bytes);

    if (!g_input_data || !g_output_data) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    mt_seed(seed);

    // Generate data with local coherence to make delta encoding meaningful
    switch (g_word_size_bytes) {
        case 1: {
            uint8_t* data = (uint8_t*)g_input_data;
            data[0] = mt_rand() & 0xFF;
            for (size_t i = 1; i < g_num_elements; ++i) {
                int16_t delta = (mt_rand() % 31) - 15; // Small delta
                data[i] = data[i - 1] + delta;
            }
            break;
        }
        case 2: {
            uint16_t* data = (uint16_t*)g_input_data;
            data[0] = mt_rand() & 0xFFFF;
            for (size_t i = 1; i < g_num_elements; ++i) {
                int32_t delta = (mt_rand() % 511) - 255;
                data[i] = data[i - 1] + delta;
            }
            break;
        }
        case 4: {
            uint32_t* data = (uint32_t*)g_input_data;
            data[0] = mt_rand();
            for (size_t i = 1; i < g_num_elements; ++i) {
                int64_t delta = (mt_rand() % 1023) - 511;
                data[i] = data[i - 1] + delta;
            }
            break;
        }
        case 8: {
            uint64_t* data = (uint64_t*)g_input_data;
            data[0] = ((uint64_t)mt_rand() << 32) | mt_rand();
            for (size_t i = 1; i < g_num_elements; ++i) {
                int64_t delta = ((int64_t)((mt_rand() % 2047) - 1023));
                data[i] = data[i - 1] + delta;
            }
            break;
        }
    }
}

void run_computation() {
    // Perform delta encoding
    switch (g_word_size_bytes) {
        case 1: {
            uint8_t* in_data = (uint8_t*)g_input_data;
            uint8_t* out_data = (uint8_t*)g_output_data;
            out_data[0] = in_data[0];
            for (size_t i = 1; i < g_num_elements; ++i) {
                out_data[i] = in_data[i] - in_data[i - 1];
            }
            break;
        }
        case 2: {
            uint16_t* in_data = (uint16_t*)g_input_data;
            uint16_t* out_data = (uint16_t*)g_output_data;
            out_data[0] = in_data[0];
            for (size_t i = 1; i < g_num_elements; ++i) {
                out_data[i] = in_data[i] - in_data[i - 1];
            }
            break;
        }
        case 4: {
            uint32_t* in_data = (uint32_t*)g_input_data;
            uint32_t* out_data = (uint32_t*)g_output_data;
            out_data[0] = in_data[0];
            for (size_t i = 1; i < g_num_elements; ++i) {
                out_data[i] = in_data[i] - in_data[i - 1];
            }
            break;
        }
        case 8: {
            uint64_t* in_data = (uint64_t*)g_input_data;
            uint64_t* out_data = (uint64_t*)g_output_data;
            out_data[0] = in_data[0];
            for (size_t i = 1; i < g_num_elements; ++i) {
                out_data[i] = in_data[i] - in_data[i - 1];
            }
            break;
        }
    }

    // Calculate a checksum to prevent dead code elimination and provide a result
    g_final_checksum = 0;
    for (size_t i = 0; i < g_total_bytes; ++i) {
        g_final_checksum += g_output_data[i];
    }
}

void cleanup() {
    free(g_input_data);
    free(g_output_data);
    g_input_data = NULL;
    g_output_data = NULL;
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%lld\n", g_final_checksum);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
