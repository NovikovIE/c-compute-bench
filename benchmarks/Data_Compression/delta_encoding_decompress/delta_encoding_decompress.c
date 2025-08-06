#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

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

// Benchmark-specific global variables
void *compressed_data = NULL;
void *decompressed_data = NULL;
size_t num_elements;
int word_size;
long long final_result = 0;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <input_size_mb> <word_size_bytes> <seed>\n", argv[0]);
        exit(1);
    }

    int input_size_mb = atoi(argv[1]);
    word_size = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (word_size != 1 && word_size != 2 && word_size != 4 && word_size != 8) {
        fprintf(stderr, "FATAL: word_size_bytes must be 1, 2, 4, or 8.\n");
        exit(1);
    }

    size_t total_bytes = (size_t)input_size_mb * 1024 * 1024;
    num_elements = total_bytes / word_size;

    compressed_data = malloc(total_bytes);
    decompressed_data = malloc(total_bytes);

    if (!compressed_data || !decompressed_data) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    mt_seed(seed);

    // Generate compressed data (a base value followed by deltas)
    switch (word_size) {
        case 1: {
            int8_t *data = (int8_t *)compressed_data;
            data[0] = (int8_t)(mt_rand() & 0xFF);
            for (size_t i = 1; i < num_elements; ++i) {
                data[i] = (int8_t)((mt_rand() % 255) - 127); // Small delta
            }
            break;
        }
        case 2: {
            int16_t *data = (int16_t *)compressed_data;
            data[0] = (int16_t)(mt_rand() & 0xFFFF);
            for (size_t i = 1; i < num_elements; ++i) {
                data[i] = (int16_t)((mt_rand() % 511) - 255); // Small delta
            }
            break;
        }
        case 4: {
            int32_t *data = (int32_t *)compressed_data;
            data[0] = (int32_t)mt_rand();
            for (size_t i = 1; i < num_elements; ++i) {
                data[i] = (int32_t)((mt_rand() % 65535) - 32767); // Small delta
            }
            break;
        }
        case 8: {
            int64_t *data = (int64_t *)compressed_data;
            data[0] = (int64_t)mt_rand() << 32 | mt_rand();
            for (size_t i = 1; i < num_elements; ++i) {
                data[i] = (int64_t)((int32_t)mt_rand()); // Moderate delta
            }
            break;
        }
    }
}

void run_computation() {
    long long checksum = 0;
    switch (word_size) {
        case 1: {
            const int8_t *c = (const int8_t *)compressed_data;
            int8_t *d = (int8_t *)decompressed_data;
            d[0] = c[0];
            checksum = d[0];
            for (size_t i = 1; i < num_elements; ++i) {
                d[i] = d[i - 1] + c[i];
                checksum += d[i];
            }
            break;
        }
        case 2: {
            const int16_t *c = (const int16_t *)compressed_data;
            int16_t *d = (int16_t *)decompressed_data;
            d[0] = c[0];
            checksum = d[0];
            for (size_t i = 1; i < num_elements; ++i) {
                d[i] = d[i - 1] + c[i];
                checksum += d[i];
            }
            break;
        }
        case 4: {
            const int32_t *c = (const int32_t *)compressed_data;
            int32_t *d = (int32_t *)decompressed_data;
            d[0] = c[0];
            checksum = d[0];
            for (size_t i = 1; i < num_elements; ++i) {
                d[i] = d[i - 1] + c[i];
                checksum += d[i];
            }
            break;
        }
        case 8: {
            const int64_t *c = (const int64_t *)compressed_data;
            int64_t *d = (int64_t *)decompressed_data;
            d[0] = c[0];
            checksum = d[0];
            for (size_t i = 1; i < num_elements; ++i) {
                d[i] = d[i - 1] + c[i];
                checksum += d[i];
            }
            break;
        }
    }
    final_result = checksum;
}

void cleanup() {
    free(compressed_data);
    free(decompressed_data);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%lld\n", final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
