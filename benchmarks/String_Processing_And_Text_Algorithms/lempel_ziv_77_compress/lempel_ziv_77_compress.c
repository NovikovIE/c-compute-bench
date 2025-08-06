#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

// --- Mersenne Twister (DO NOT MODIFY) ---
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
// --- End Mersenne Twister ---

// Benchmark-specific constants
#define VOCAB_SIZE 256
#define MAX_WORD_LEN 16
#define MIN_WORD_LEN 4
#define MIN_MATCH_LEN 3
#define MAX_LOOKAHEAD 258 // As in DEFLATE

// Global data structure
struct {
    char* input_data;
    size_t input_data_size;
    int sliding_window_size;
    
    // Output is a series of (distance, length, next_char) tuples
    int* compressed_output;
    size_t compressed_token_count;

    // Final result to prevent dead code elimination
    long long result_checksum;
} g_data;

// Generate input data with repeating patterns to make compression meaningful
void generate_repetitive_text(char* buffer, size_t size) {
    char vocabulary[VOCAB_SIZE][MAX_WORD_LEN + 1];
    int word_lengths[VOCAB_SIZE];

    // Create a vocabulary of random words
    for (int i = 0; i < VOCAB_SIZE; i++) {
        word_lengths[i] = (mt_rand() % (MAX_WORD_LEN - MIN_WORD_LEN + 1)) + MIN_WORD_LEN;
        for (int j = 0; j < word_lengths[i]; j++) {
            vocabulary[i][j] = (char)((mt_rand() % 94) + 33); // Printable ASCII
        }
        vocabulary[i][word_lengths[i]] = '\0';
    }

    // Fill the buffer by randomly picking words from the vocabulary
    size_t current_pos = 0;
    while (current_pos < size) {
        int word_idx = mt_rand() % VOCAB_SIZE;
        int len_to_copy = word_lengths[word_idx];
        if (current_pos + len_to_copy >= size) {
            len_to_copy = size - current_pos;
        }
        if (len_to_copy > 0) {
            memcpy(buffer + current_pos, vocabulary[word_idx], len_to_copy);
        }
        current_pos += len_to_copy;
    }
}


void setup_benchmark(int argc, char* argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <input_data_size_kb> <sliding_window_size> <seed>\n", argv[0]);
        exit(1);
    }

    int input_data_size_kb = atoi(argv[1]);
    g_data.sliding_window_size = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    mt_seed(seed);
    
    g_data.input_data_size = (size_t)input_data_size_kb * 1024;
    if (g_data.sliding_window_size <= 0 || g_data.input_data_size <= 0) {
        fprintf(stderr, "FATAL: Invalid size parameters.\n");
        exit(1);
    }

    g_data.input_data = (char*)malloc(g_data.input_data_size);
    if (!g_data.input_data) {
        fprintf(stderr, "FATAL: Memory allocation failed for input_data.\n");
        exit(1);
    }
    
    // Worst case: no matches, one (0,0,char) token per byte. A token is 3 ints.
    size_t max_tokens = g_data.input_data_size;
    g_data.compressed_output = (int*)malloc(max_tokens * 3 * sizeof(int));
     if (!g_data.compressed_output) {
        fprintf(stderr, "FATAL: Memory allocation failed for compressed_output.\n");
        free(g_data.input_data);
        exit(1);
    }

    generate_repetitive_text(g_data.input_data, g_data.input_data_size);
    
    g_data.result_checksum = 0;
    g_data.compressed_token_count = 0;
}

void run_computation() {
    size_t cursor = 0;
    size_t output_idx = 0;
    g_data.result_checksum = 0;

    while (cursor < g_data.input_data_size) {
        int best_len = 0;
        int best_dist = 0;

        size_t search_start = (cursor > (size_t)g_data.sliding_window_size) ? (cursor - g_data.sliding_window_size) : 0;
        
        // Search backwards from current position for the longest match
        for (size_t search_pos = cursor - 1; search_pos != (size_t)-1 && search_pos >= search_start; --search_pos) {
            int current_len = 0;
            // Compare characters to find match length
            while (cursor + current_len < g_data.input_data_size &&
                   g_data.input_data[search_pos + current_len] == g_data.input_data[cursor + current_len] &&
                   current_len < MAX_LOOKAHEAD) {
                current_len++;
            }

            if (current_len > best_len) {
                best_len = current_len;
                best_dist = cursor - search_pos;
            }
        }

        if (best_len >= MIN_MATCH_LEN) {
            // Found a match. The "next char" is the char after the match in the lookahead buffer.
            char next_char = (cursor + best_len < g_data.input_data_size) 
                             ? g_data.input_data[cursor + best_len] 
                             : '\0';
            
            g_data.compressed_output[output_idx++] = best_dist;
            g_data.compressed_output[output_idx++] = best_len;
            g_data.compressed_output[output_idx++] = (int)next_char;

            g_data.result_checksum += best_dist + best_len + (int)next_char;
            cursor += best_len + 1; // Advance cursor past the match and the next char
        } else {
            // No match found, emit a literal
            char literal = g_data.input_data[cursor];
            
            g_data.compressed_output[output_idx++] = 0; // distance 0 indicates literal
            g_data.compressed_output[output_idx++] = 0; // length 0 indicates literal
            g_data.compressed_output[output_idx++] = (int)literal;

            g_data.result_checksum += (int)literal;
            cursor++; // Advance cursor by one
        }
    }
    g_data.compressed_token_count = output_idx / 3;
}


void cleanup() {
    free(g_data.input_data);
    free(g_data.compressed_output);
}

int main(int argc, char* argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);
    
    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
    // Print result to stdout
    printf("%lld\n", g_data.result_checksum);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}
