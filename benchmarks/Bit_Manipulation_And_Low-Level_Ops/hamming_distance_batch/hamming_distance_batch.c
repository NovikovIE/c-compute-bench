#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <inttypes.h> // For PRIu64

// --- Mersenne Twister (MT19937) Generator ---
// Provided verbatim as per instructions
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

// --- Benchmark Globals ---

// Parameters
int NUM_PAIRS;
int BIT_LENGTH;

// Derived values
size_t NUM_WORDS_PER_SEQUENCE;
size_t TOTAL_WORDS;

// Data arrays
uint64_t *data_a = NULL;
uint64_t *data_b = NULL;

// Result accumulator
uint64_t total_hamming_distance = 0;


// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_pairs> <bit_length> <seed>\n", argv[0]);
        exit(1);
    }

    NUM_PAIRS = atoi(argv[1]);
    BIT_LENGTH = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);
    
    if (NUM_PAIRS <= 0 || BIT_LENGTH <= 0) {
        fprintf(stderr, "ERROR: num_pairs and bit_length must be positive integers.\n");
        exit(1);
    }

    // Calculate memory requirements
    NUM_WORDS_PER_SEQUENCE = (BIT_LENGTH + 63) / 64;
    TOTAL_WORDS = (size_t)NUM_PAIRS * NUM_WORDS_PER_SEQUENCE;

    // Allocate memory
    data_a = (uint64_t*)malloc(TOTAL_WORDS * sizeof(uint64_t));
    if (!data_a) {
        fprintf(stderr, "FATAL: Memory allocation for data_a failed.\n");
        exit(1);
    }

    data_b = (uint64_t*)malloc(TOTAL_WORDS * sizeof(uint64_t));
    if (!data_b) {
        fprintf(stderr, "FATAL: Memory allocation for data_b failed.\n");
        free(data_a);
        exit(1);
    }
    
    // Seed the random number generator
    mt_seed(seed);

    // Generate random bit sequences
    // We combine two 32-bit random numbers to create one 64-bit number
    for (size_t i = 0; i < TOTAL_WORDS; ++i) {
        data_a[i] = ((uint64_t)mt_rand() << 32) | mt_rand();
        data_b[i] = ((uint64_t)mt_rand() << 32) | mt_rand();
    }
    
    // If bit_length is not a multiple of 64, we must mask the last word of each sequence
    // to ensure we only consider the specified number of bits.
    int last_word_bits = BIT_LENGTH % 64;
    if (last_word_bits > 0) {
        uint64_t mask = (1ULL << last_word_bits) - 1;
        for (int i = 0; i < NUM_PAIRS; ++i) {
            size_t last_word_index = (i + 1) * NUM_WORDS_PER_SEQUENCE - 1;
            data_a[last_word_index] &= mask;
            data_b[last_word_index] &= mask;
        }
    }
}

void run_computation() {
    total_hamming_distance = 0; // Reset before computation
    for (int i = 0; i < NUM_PAIRS; ++i) {
        size_t start_index = i * NUM_WORDS_PER_SEQUENCE;
        for (size_t j = 0; j < NUM_WORDS_PER_SEQUENCE; ++j) {
            // XOR the two words to find the differing bits
            uint64_t xor_result = data_a[start_index + j] ^ data_b[start_index + j];
            
            // Count the number of set bits (1s) in the XOR result.
            // This is the Hamming distance for this 64-bit chunk.
            // __builtin_popcountll is a highly optimized GCC/Clang intrinsic
            // that often compiles to a single CPU instruction (e.g., popcnt).
            total_hamming_distance += __builtin_popcountll(xor_result);
        }
    }
}

void cleanup() {
    free(data_a);
    free(data_b);
    data_a = NULL;
    data_b = NULL;
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
    // Print the final accumulated Hamming distance to stdout
    printf("%" PRIu64 "\n", total_hamming_distance);

    // Print the execution time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
