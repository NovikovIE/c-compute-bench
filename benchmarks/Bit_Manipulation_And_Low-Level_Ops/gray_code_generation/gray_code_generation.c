#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <inttypes.h>

// --- Mersenne Twister (MT19937) Generator ---
// Do Not Modify
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
// --- End of Mersenne Twister ---

#define DUMMY_ALLOC_SIZE 1024

// --- Benchmark Globals ---
int num_bits;
uint64_t count;
uint64_t accumulated_sum;
uint32_t* dummy_data; // Unused data to demonstrate setup/cleanup memory management

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_bits> <seed>\n", argv[0]);
        exit(1);
    }

    num_bits = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if (num_bits < 1 || num_bits > 32) {
        fprintf(stderr, "Error: num_bits must be between 1 and 32.\n");
        exit(1);
    }

    // Seed the random number generator
    mt_seed(seed);
    
    // Calculate the number of iterations for the computation
    count = 1ULL << num_bits;
    accumulated_sum = 0;

    // Allocate and initialize some dummy data to adhere to the heap allocation requirement.
    // This is part of setup and not timed.
    dummy_data = (uint32_t*)malloc(DUMMY_ALLOC_SIZE * sizeof(uint32_t));
    if (dummy_data == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }
    for (int i = 0; i < DUMMY_ALLOC_SIZE; ++i) {
        dummy_data[i] = mt_rand();
    }
}

void run_computation() {
    uint64_t sum = 0;
    // The loop iterates 2^num_bits times.
    for (uint64_t i = 0; i < count; ++i) {
         // The core operation: generate the Gray code for 'i' and add it to the sum.
         // The binary-reflected Gray code for an integer 'n' is 'n' XOR ('n' shifted right by 1).
        sum += (i ^ (i >> 1));
    }
    accumulated_sum = sum;
}

void cleanup() {
    // Free any memory allocated in setup_benchmark
    free(dummy_data);
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

    // Print the final accumulated result to stdout to prevent dead code elimination
    printf("%" PRIu64 "\n", accumulated_sum);

    // Print the execution time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
