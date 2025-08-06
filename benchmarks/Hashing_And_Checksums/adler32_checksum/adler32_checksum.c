#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// Mersenne Twister (verbatim)
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
// End Mersenne Twister

// Benchmark Globals
unsigned char* data_buffer;
size_t data_size;
uint32_t final_checksum;

#define ADLER32_MOD 65521

// Setup: parse args, allocate and fill buffer
void setup_benchmark(int argc, char* argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <data_size_kb> <seed>\n", argv[0]);
        exit(1);
    }
    
    long data_size_kb = atol(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    data_size = (size_t)data_size_kb * 1024;
    
    mt_seed(seed);
    
    data_buffer = (unsigned char*)malloc(data_size);
    if (data_buffer == NULL) {
        fprintf(stderr, "Failed to allocate memory for data buffer.\n");
        exit(1);
    }
    
    for (size_t i = 0; i < data_size; i++) {
        data_buffer[i] = mt_rand() & 0xFF;
    }
}

// Computation: calculate Adler-32 checksum
void run_computation() {
    uint32_t s1 = 1;
    uint32_t s2 = 0;

    for (size_t i = 0; i < data_size; ++i) {
        s1 = (s1 + data_buffer[i]) % ADLER32_MOD;
        s2 = (s2 + s1) % ADLER32_MOD;
    }

    final_checksum = (s2 << 16) | s1;
}

// Cleanup: free allocated memory
void cleanup() {
    if (data_buffer != NULL) {
        free(data_buffer);
    }
}


// Main function: orchestrate and time the benchmark
int main(int argc, char* argv[]) {
    struct timespec start, end;
    
    setup_benchmark(argc, argv);
    
    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);
    
    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
    // Print the final checksum to stdout
    printf("%u\n", final_checksum);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);
    
    return 0;
}
