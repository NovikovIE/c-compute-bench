#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- Mersenne Twister (MT19937) --- (DO NOT MODIFY)
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
// --- End of Mersenne Twister ---

// --- Benchmark Globals ---
uint64_t* source_matrix = NULL;
uint64_t* dest_matrix = NULL;
unsigned int matrix_size;
unsigned long long final_checksum = 0;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <matrix_size> <seed>\n", argv[0]);
        exit(1);
    }

    matrix_size = atoi(argv[1]);
    uint32_t seed = atoi(argv[2]);

    if (matrix_size <= 0) {
        fprintf(stderr, "FATAL: matrix_size must be a positive integer.\n");
        exit(1);
    }

    mt_seed(seed);

    uint64_t num_bits = (uint64_t)matrix_size * matrix_size;
    size_t num_words = (num_bits + 63) / 64;

    source_matrix = (uint64_t*)malloc(num_words * sizeof(uint64_t));
    if (!source_matrix) {
        fprintf(stderr, "FATAL: Failed to allocate memory for source matrix.\n");
        exit(1);
    }

    // Use calloc to zero-initialize the destination matrix
    dest_matrix = (uint64_t*)calloc(num_words, sizeof(uint64_t));
    if (!dest_matrix) {
        fprintf(stderr, "FATAL: Failed to allocate memory for destination matrix.\n");
        free(source_matrix);
        exit(1);
    }

    // Populate source matrix with random bits
    for (size_t i = 0; i < num_words; ++i) {
        source_matrix[i] = ((uint64_t)mt_rand() << 32) | mt_rand();
    }
}

void run_computation() {
    // Naive bit-by-bit matrix transposition
    for (unsigned int r = 0; r < matrix_size; ++r) {
        for (unsigned int c = 0; c < matrix_size; ++c) {
            // Calculate source bit position
            uint64_t src_bit_idx = (uint64_t)r * matrix_size + c;
            uint64_t src_word_idx = src_bit_idx >> 6;      // same as / 64
            uint64_t src_bit_in_word = src_bit_idx & 63; // same as % 64

            // Extract the bit
            uint64_t bit = (source_matrix[src_word_idx] >> src_bit_in_word) & 1ULL;

            // If the bit is 1, set the corresponding bit in the destination
            if (bit) {
                // Calculate destination bit position (transposed)
                uint64_t dest_bit_idx = (uint64_t)c * matrix_size + r;
                uint64_t dest_word_idx = dest_bit_idx >> 6;
                uint64_t dest_bit_in_word = dest_bit_idx & 63;
                dest_matrix[dest_word_idx] |= (1ULL << dest_bit_in_word);
            }
        }
    }

    // Calculate a checksum of the transposed matrix to prevent dead code elimination
    uint64_t num_bits = (uint64_t)matrix_size * matrix_size;
    size_t num_words = (num_bits + 63) / 64;
    unsigned long long checksum = 0;
    for (size_t i = 0; i < num_words; ++i) {
        checksum ^= dest_matrix[i];
    }
    final_checksum = checksum;
}

void cleanup() {
    if (source_matrix) {
        free(source_matrix);
        source_matrix = NULL;
    }
    if (dest_matrix) {
        free(dest_matrix);
        dest_matrix = NULL;
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

    // Print result to stdout
    printf("%llu\n", final_checksum);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
