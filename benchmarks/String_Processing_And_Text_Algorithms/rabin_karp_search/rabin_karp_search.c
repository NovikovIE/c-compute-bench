#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

// --- BEGIN: Mersenne Twister (DO NOT MODIFY) ---
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
// --- END: Mersenne Twister ---

// --- Benchmark Globals ---
static size_t text_length;
static size_t pattern_length;
static char* text;
static char* pattern;
static int total_matches; // Accumulated result

// Rabin-Karp constants
#define ALPHABET_SIZE 4 // Using a small alphabet 'A', 'B', 'C', 'D'
#define Q 1000000007  // A large prime number for modulo operations

// --- Benchmark Functions ---

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <text_length> <pattern_length> <seed>\n", argv[0]);
        exit(1);
    }

    text_length = atol(argv[1]);
    pattern_length = atol(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (pattern_length > text_length) {
        fprintf(stderr, "Error: pattern_length cannot be greater than text_length.\n");
        exit(1);
    }

    mt_seed(seed);

    text = (char*)malloc(text_length + 1);
    pattern = (char*)malloc(pattern_length + 1);

    if (!text || !pattern) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(1);
    }

    for (size_t i = 0; i < text_length; i++) {
        text[i] = 'A' + (mt_rand() % ALPHABET_SIZE);
    }
    text[text_length] = '\0';

    for (size_t i = 0; i < pattern_length; i++) {
        pattern[i] = 'A' + (mt_rand() % ALPHABET_SIZE);
    }
    pattern[pattern_length] = '\0';

    total_matches = 0;
}

void run_computation() {
    long long p_hash = 0; // hash value for pattern
    long long t_hash = 0; // hash value for text window
    long long h = 1;      // (ALPHABET_SIZE^(pattern_length-1)) % Q

    // Calculate h = ALPHABET_SIZE^(pattern_length-1) % Q
    for (size_t i = 0; i < pattern_length - 1; i++) {
        h = (h * ALPHABET_SIZE) % Q;
    }

    // Calculate the hash value of pattern and first window of text
    for (size_t i = 0; i < pattern_length; i++) {
        p_hash = (ALPHABET_SIZE * p_hash + (pattern[i] - 'A')) % Q;
        t_hash = (ALPHABET_SIZE * t_hash + (text[i] - 'A')) % Q;
    }

    // Slide the pattern over text one by one
    size_t limit = text_length - pattern_length;
    for (size_t i = 0; i <= limit; i++) {
        // Check the hash values of current window of text and pattern
        // If the hash values match then only check for characters one by one
        if (p_hash == t_hash) {
            // Check for characters one by one (spurious hit check)
            if (memcmp(pattern, text + i, pattern_length) == 0) {
                total_matches++;
            }
        }

        // Calculate hash value for next window of text: 
        // Remove leading digit, add trailing digit
        if (i < limit) {
            long long leading_char_val = (long long)(text[i] - 'A') * h;
            t_hash = (ALPHABET_SIZE * (t_hash - leading_char_val) + (text[i + pattern_length] - 'A')) % Q;

            // We might get a negative hash value, convert it to positive
            if (t_hash < 0) {
                t_hash = (t_hash + Q);
            }
        }
    }
}

void cleanup() {
    free(text);
    free(pattern);
}

int main(int argc, char* argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%d\n", total_matches);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
