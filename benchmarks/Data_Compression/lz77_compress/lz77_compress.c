#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- Mersenne Twister (MT19937) Vebatim --- 
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

// Type definition for an LZ77 compression token
typedef struct {
    int offset; // Offset back into the search buffer
    int length; // Length of the match
    unsigned char literal; // The literal character if no match is found
} lz77_token_t;

// Global data structure to hold benchmark data and parameters
typedef struct {
    unsigned char *input_data;
    size_t input_size;

    lz77_token_t *compressed_tokens;
    size_t max_tokens;
    size_t token_count; // Final result: number of tokens generated

    // Parameters from command line
    size_t sliding_window_size;
    size_t lookahead_buffer_size;
} benchmark_data_t;

static benchmark_data_t g_data;
static const int MIN_MATCH_LENGTH = 3;

// Setup: parse arguments, allocate memory, generate input data
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <input_size_mb> <sliding_window_size> <lookahead_buffer_size> <seed>\n", argv[0]);
        exit(1);
    }

    long input_size_mb = atol(argv[1]);
    g_data.sliding_window_size = atol(argv[2]);
    g_data.lookahead_buffer_size = atol(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    if (input_size_mb <= 0 || g_data.sliding_window_size <= 0 || g_data.lookahead_buffer_size <= 0) {
        fprintf(stderr, "FATAL: All size parameters must be positive.\n");
        exit(1);
    }

    g_data.input_size = input_size_mb * 1024 * 1024;

    mt_seed(seed);

    // Allocate input data buffer
    g_data.input_data = (unsigned char *)malloc(g_data.input_size);
    if (!g_data.input_data) {
        fprintf(stderr, "FATAL: Failed to allocate memory for input data\n");
        exit(1);
    }

    // Generate random input data
    size_t i = 0;
    while (i < g_data.input_size) {
        uint32_t r = mt_rand();
        g_data.input_data[i++] = (unsigned char)(r & 0xFF);
        if (i < g_data.input_size) g_data.input_data[i++] = (unsigned char)((r >> 8) & 0xFF);
        if (i < g_data.input_size) g_data.input_data[i++] = (unsigned char)((r >> 16) & 0xFF);
        if (i < g_data.input_size) g_data.input_data[i++] = (unsigned char)((r >> 24) & 0xFF);
    }

    // Allocate space for compressed tokens. Worst case is one token per byte.
    g_data.max_tokens = g_data.input_size;
    g_data.compressed_tokens = (lz77_token_t *)malloc(g_data.max_tokens * sizeof(lz77_token_t));
    if (!g_data.compressed_tokens) {
        fprintf(stderr, "FATAL: Failed to allocate memory for compressed tokens\n");
        free(g_data.input_data);
        exit(1);
    }

    g_data.token_count = 0;
}

// Computation: perform LZ77 compression
void run_computation() {
    size_t cursor = 0;

    while (cursor < g_data.input_size) {
        int best_len = 0;
        int best_off = 0;

        size_t search_start = (cursor > g_data.sliding_window_size) ? (cursor - g_data.sliding_window_size) : 0;

        // Search for the longest match in the sliding window
        // The loop is safe for size_t wraparound because i will never be < search_start before it wraps.
        for (size_t i = cursor - 1; i < cursor && i >= search_start; i--) {
            int current_len = 0;
            while (cursor + current_len < g_data.input_size &&
                   g_data.input_data[i + current_len] == g_data.input_data[cursor + current_len] &&
                   current_len < g_data.lookahead_buffer_size) {
                current_len++;
            }

            if (current_len > best_len) {
                best_len = current_len;
                best_off = cursor - i;
                if (best_len == g_data.lookahead_buffer_size) {
                    break; // Max possible match found, exit search early
                }
            }
        }

        if (g_data.token_count >= g_data.max_tokens) {
            // Should not happen with current allocation strategy, but good practice
            break;
        }

        // If a sufficiently long match is found, create a pointer token
        if (best_len >= MIN_MATCH_LENGTH) {
            g_data.compressed_tokens[g_data.token_count].offset = best_off;
            g_data.compressed_tokens[g_data.token_count].length = best_len;
            g_data.compressed_tokens[g_data.token_count].literal = 0;
            g_data.token_count++;
            cursor += best_len;
        } else {
            // Otherwise, create a literal token
            g_data.compressed_tokens[g_data.token_count].offset = 0;
            g_data.compressed_tokens[g_data.token_count].length = 0;
            g_data.compressed_tokens[g_data.token_count].literal = g_data.input_data[cursor];
            g_data.token_count++;
            cursor++;
        }
    }
}

// Cleanup: free all allocated memory
void cleanup() {
    free(g_data.input_data);
    g_data.input_data = NULL;
    free(g_data.compressed_tokens);
    g_data.compressed_tokens = NULL;
}

int main(int argc, char *argv[]) {
    struct timespec start, end;
    size_t final_result;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    final_result = g_data.token_count;

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final result (number of tokens) to stdout
    printf("%zu\n", final_result);

    // Print the execution time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}