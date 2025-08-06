#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <string.h>

// Mersenne Twister (MT19937) Generator - DO NOT MODIFY
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
// End of Mersenne Twister Generator

// --- Benchmark Globals ---
size_t sequence_length;
int max_pattern_size;
char* dna_sequence;
unsigned long long total_repeats_found;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <sequence_length> <max_pattern_size> <seed>\n", argv[0]);
        exit(1);
    }

    sequence_length = (size_t)atol(argv[1]);
    max_pattern_size = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    mt_seed(seed);

    dna_sequence = (char*)malloc((sequence_length + 1) * sizeof(char));
    if (dna_sequence == NULL) {
        fprintf(stderr, "Memory allocation failed for DNA sequence.\n");
        exit(1);
    }

    const char bases[] = "ACGT";
    for (size_t i = 0; i < sequence_length; ++i) {
        dna_sequence[i] = bases[mt_rand() % 4];
    }
    dna_sequence[sequence_length] = '\0';
}

void run_computation() {
    total_repeats_found = 0;

    // Iterate through all possible start positions in the sequence
    for (size_t i = 0; i < sequence_length; ++i) {
        // Iterate through all possible pattern lengths from 2 up to max_pattern_size
        for (int p_len = 2; p_len <= max_pattern_size; ++p_len) {
            // Check if there is enough sequence left for a pattern and its adjacent repeat
            if (i + 2 * p_len > sequence_length) {
                break; // Not enough space for this pattern length onwards
            }
            
            // Compare the block of size p_len at 'i' with the adjacent block at 'i + p_len'
            if (strncmp(&dna_sequence[i], &dna_sequence[i + p_len], p_len) == 0) {
                total_repeats_found++;
            }
        }
    }
}

void cleanup() {
    free(dna_sequence);
}

// --- Main and timing ---

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout
    printf("%llu\n", total_repeats_found);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
