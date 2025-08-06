#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// Embedded Mersenne Twister Begin
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
// Embedded Mersenne Twister End

// Benchmark parameters
static size_t text_length;
static size_t pattern_length;
static int alphabet_size;

// Data structures
static char *text = NULL;
static char *pattern = NULL;
static int *bad_char_table = NULL;

// Result
static int total_matches = 0;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s text_length pattern_length alphabet_size seed\n", argv[0]);
        exit(1);
    }

    text_length = atol(argv[1]);
    pattern_length = atol(argv[2]);
    alphabet_size = atoi(argv[3]);
    uint32_t seed = atoi(argv[4]);

    if (pattern_length == 0 || text_length == 0) {
        fprintf(stderr, "FATAL: Text and pattern lengths must be positive.\n");
        exit(1);
    }
    if (pattern_length > text_length) {
        fprintf(stderr, "FATAL: Pattern length cannot be greater than text length.\n");
        exit(1);
    }
    if (alphabet_size <= 0 || alphabet_size > 256) {
        fprintf(stderr, "FATAL: Alphabet size must be between 1 and 256.\n");
        exit(1);
    }

    mt_seed(seed);

    text = (char *)malloc(text_length * sizeof(char));
    pattern = (char *)malloc(pattern_length * sizeof(char));
    bad_char_table = (int *)malloc(alphabet_size * sizeof(int));

    if (!text || !pattern || !bad_char_table) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    for (size_t i = 0; i < text_length; ++i) {
        text[i] = (char)(mt_rand() % alphabet_size);
    }

    for (size_t i = 0; i < pattern_length; ++i) {
        pattern[i] = (char)(mt_rand() % alphabet_size);
    }
}

void run_computation() {
    // Preprocessing phase: Create the bad character table
    for (int i = 0; i < alphabet_size; i++) {
        bad_char_table[i] = -1;
    }
    for (size_t i = 0; i < pattern_length; i++) {
        bad_char_table[(unsigned char)pattern[i]] = (int)i;
    }

    // Search phase
    size_t s = 0; // s is the shift of the pattern with respect to the text
    while (s <= (text_length - pattern_length)) {
        int j = pattern_length - 1;

        while (j >= 0 && pattern[j] == text[s + j]) {
            j--;
        }

        if (j < 0) {
            total_matches++;
            // To find all occurrences, we must shift by at least 1.
            // A more advanced shift can be used here, but for simplicity
            // and to guarantee progress, we use a simple shift after a match.
            s++;
        } else {
            // Shift based on the bad character rule.
            // The cast to unsigned char is crucial to prevent negative indexing.
            int bad_char_pos = bad_char_table[(unsigned char)text[s + j]];
            int shift = j - bad_char_pos;
            // Ensure we make positive progress.
            s += (shift > 1) ? shift : 1;
        }
    }
}

void cleanup() {
    free(text);
    free(pattern);
    free(bad_char_table);
    text = NULL;
    pattern = NULL;
    bad_char_table = NULL;
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%d\n", total_matches);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
