/*
 * lzss_compress: A benchmark for LZSS data compression.
 * 
 * This program generates semi-compressible random data and then compresses it
 * using a simplified Lempel-Ziv-Storer-Szymanski (LZSS) algorithm. The core
 * of the computation is the search for repeating byte sequences within a
 * sliding window, which is CPU-intensive.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

// --- BEGIN MERSENNE TWISTER (MT19937) --- Do Not Modify ---
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
// --- END MERSENNE TWISTER ---

#define MIN_MATCH_LENGTH 3
#define MAX_MATCH_LENGTH 255

// Use a global struct to hold benchmark data
typedef struct {
    unsigned char *input_data;
    unsigned char *compressed_data;
    size_t input_size;
    size_t compressed_capacity;
    int window_size;
} BenchmarkData;

BenchmarkData g_data;
size_t final_result; // Accumulated result to prevent dead code elimination

// Data setup: parse arguments, allocate memory, generate input data
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <input_size_mb> <sliding_window_size> <seed>\n", argv[0]);
        exit(1);
    }

    long input_size_mb = atol(argv[1]);
    g_data.window_size = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    mt_seed(seed);

    g_data.input_size = input_size_mb * 1024 * 1024;

    if (g_data.window_size <= 0 || g_data.window_size > 65535) {
        fprintf(stderr, "FATAL: sliding_window_size must be between 1 and 65535\n");
        exit(1);
    }

    g_data.input_data = (unsigned char *)malloc(g_data.input_size);
    if (!g_data.input_data) {
        fprintf(stderr, "FATAL: Failed to allocate memory for input data\n");
        exit(1);
    }

    // Compressed data can be larger in worst case, allocate generously
    g_data.compressed_capacity = g_data.input_size * 2;
    g_data.compressed_data = (unsigned char *)malloc(g_data.compressed_capacity);
    if (!g_data.compressed_data) {
        fprintf(stderr, "FATAL: Failed to allocate memory for compressed data\n");
        exit(1);
    }

    // Generate semi-compressible data by repeating random "words"
    #define DICTIONARY_WORDS 512
    #define MAX_WORD_LEN 128
    #define MIN_WORD_LEN 4
    unsigned char* dictionary[DICTIONARY_WORDS];
    int word_lengths[DICTIONARY_WORDS];

    for (int i = 0; i < DICTIONARY_WORDS; ++i) {
        word_lengths[i] = (mt_rand() % (MAX_WORD_LEN - MIN_WORD_LEN + 1)) + MIN_WORD_LEN;
        dictionary[i] = (unsigned char*) malloc(word_lengths[i]);
        for (int j = 0; j < word_lengths[i]; ++j) {
            dictionary[i][j] = mt_rand() % 256;
        }
    }

    size_t current_pos = 0;
    while (current_pos < g_data.input_size) {
        int word_idx = mt_rand() % DICTIONARY_WORDS;
        int len_to_copy = word_lengths[word_idx];
        if (current_pos + len_to_copy > g_data.input_size) {
            len_to_copy = g_data.input_size - current_pos;
        }
        memcpy(g_data.input_data + current_pos, dictionary[word_idx], len_to_copy);
        current_pos += len_to_copy;
    }

    for (int i = 0; i < DICTIONARY_WORDS; ++i) {
        free(dictionary[i]);
    }
}

// Core computation: LZSS compression algorithm
void run_computation() {
    size_t input_pos = 0;
    size_t output_pos = 0;

    while (input_pos < g_data.input_size) {
        int best_match_len = 0;
        int best_match_offset = 0;

        size_t search_start = (input_pos > g_data.window_size) ? (input_pos - g_data.window_size) : 0;

        for (size_t i = search_start; i < input_pos; ++i) {
            int current_len = 0;
            int max_possible_len = g_data.input_size - input_pos;
            if (max_possible_len > MAX_MATCH_LENGTH) {
                max_possible_len = MAX_MATCH_LENGTH;
            }
            while (current_len < max_possible_len && g_data.input_data[i + current_len] == g_data.input_data[input_pos + current_len]) {
                current_len++;
            }

            if (current_len > best_match_len) {
                best_match_len = current_len;
                best_match_offset = input_pos - i;
            }
        }

        if (output_pos + 4 >= g_data.compressed_capacity) {
            // Prevent buffer overflow, though our capacity should be sufficient
            break;
        }

        if (best_match_len >= MIN_MATCH_LENGTH) {
            // Write a simple token: <length, offset>
            g_data.compressed_data[output_pos++] = (unsigned char)best_match_len;
            g_data.compressed_data[output_pos++] = (unsigned char)(best_match_offset >> 8);
            g_data.compressed_data[output_pos++] = (unsigned char)(best_match_offset & 0xFF);
            input_pos += best_match_len;
        } else {
            // Write a literal byte, marked by a 0 length token
            g_data.compressed_data[output_pos++] = 0;
            g_data.compressed_data[output_pos++] = g_data.input_data[input_pos];
            input_pos++;
        }
    }

    final_result = output_pos;
}

// Free all allocated memory
void cleanup() {
    free(g_data.input_data);
    free(g_data.compressed_data);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final result (compressed size) to stdout
    printf("%zu\n", final_result);

    // Print the execution time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
