#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- Mersenne Twister (Verbatim) ---
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

// --- Global Data ---
unsigned char* data_buffer = NULL;
size_t data_buffer_size = 0;
uint32_t final_hash = 0;

// --- Benchmark Algorithm ---
uint32_t jenkins_one_at_a_time_hash(const uint8_t* key, size_t length) {
    size_t i = 0;
    uint32_t hash = 0;
    while (i != length) {
        hash += key[i++];
        hash += hash << 10;
        hash ^= hash >> 6;
    }
    hash += hash << 3;
    hash ^= hash >> 11;
    hash += hash << 15;
    return hash;
}

// --- Benchmark Functions ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s data_size_kb seed\n", argv[0]);
        exit(1);
    }

    long data_size_kb = atol(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    data_buffer_size = data_size_kb * 1024;

    mt_seed(seed);

    data_buffer = (unsigned char *)malloc(data_buffer_size);
    if (data_buffer == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Populate the buffer with random 32-bit integers for efficiency
    uint32_t* p = (uint32_t*)data_buffer;
    size_t num_u32 = data_buffer_size / sizeof(uint32_t);
    for (size_t i = 0; i < num_u32; ++i) {
        p[i] = mt_rand();
    }

    // Fill remaining bytes if size is not a multiple of 4
    size_t remainder_start = num_u32 * sizeof(uint32_t);
    for (size_t i = remainder_start; i < data_buffer_size; ++i) {
        data_buffer[i] = (unsigned char)(mt_rand() & 0xff);
    }
}

void run_computation() {
    final_hash = jenkins_one_at_a_time_hash(data_buffer, data_buffer_size);
}

void cleanup() {
    if (data_buffer != NULL) {
        free(data_buffer);
        data_buffer = NULL;
    }
}

// --- Main Execution ---
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%u\n", final_hash);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
