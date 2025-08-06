#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- Mersenne Twister (MT19937) Generator ---
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

// --- Benchmark Data and Globals ---
typedef struct {
    size_t text_length;
    size_t pattern_length;
    char *text;
    char *pattern;
    int *lps; // Longest Proper Prefix Suffix array
    int match_count; // Result of the computation
} BenchmarkData;

static BenchmarkData g_data;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s text_length pattern_length seed\n", argv[0]);
        exit(1);
    }

    g_data.text_length = atol(argv[1]);
    g_data.pattern_length = atol(argv[2]);
    uint32_t seed = atoi(argv[3]);

    if (g_data.pattern_length == 0 || g_data.text_length == 0 || g_data.pattern_length > g_data.text_length) {
        fprintf(stderr, "FATAL: Invalid lengths. Must have text_length > 0, pattern_length > 0, and text_length >= pattern_length.\n");
        exit(1);
    }

    mt_seed(seed);

    g_data.text = (char*)malloc(g_data.text_length * sizeof(char));
    g_data.pattern = (char*)malloc(g_data.pattern_length * sizeof(char));
    g_data.lps = (int*)malloc(g_data.pattern_length * sizeof(int));

    if (!g_data.text || !g_data.pattern || !g_data.lps) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Use a small alphabet to increase partial matches, exercising the KMP failure function.
    const char alphabet[] = "ABCD";
    const int alphabet_size = sizeof(alphabet) - 1;

    for (size_t i = 0; i < g_data.text_length; i++) {
        g_data.text[i] = alphabet[mt_rand() % alphabet_size];
    }

    for (size_t i = 0; i < g_data.pattern_length; i++) {
        g_data.pattern[i] = alphabet[mt_rand() % alphabet_size];
    }
}

// Computes the Longest Proper Prefix which is also a Suffix (LPS) array.
void compute_lps_array() {
    size_t length = 0;
    g_data.lps[0] = 0;
    size_t i = 1;

    while (i < g_data.pattern_length) {
        if (g_data.pattern[i] == g_data.pattern[length]) {
            length++;
            g_data.lps[i] = length;
            i++;
        } else {
            if (length != 0) {
                length = g_data.lps[length - 1];
            } else {
                g_data.lps[i] = 0;
                i++;
            }
        }
    }
}

void run_computation() {
    // Phase 1: Preprocess the pattern to build the LPS array.
    compute_lps_array();

    // Phase 2: Search for the pattern in the text using the LPS array.
    size_t i = 0; // index for text
    size_t j = 0; // index for pattern
    int count = 0;

    while (i < g_data.text_length) {
        if (g_data.pattern[j] == g_data.text[i]) {
            i++;
            j++;
        }

        if (j == g_data.pattern_length) {
            count++;
            // Found a match, continue searching for next match.
            j = g_data.lps[j - 1];
        } else if (i < g_data.text_length && g_data.pattern[j] != g_data.text[i]) {
            // Mismatch after a partial match. Use LPS array to avoid redundant comparisons.
            if (j != 0) {
                j = g_data.lps[j - 1];
            } else {
                // Mismatch at the beginning of the pattern, just advance text pointer.
                i++;
            }
        }
    }
    g_data.match_count = count;
}

void cleanup() {
    free(g_data.text);
    free(g_data.pattern);
    free(g_data.lps);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print accumulated result to stdout to prevent dead code elimination
    printf("%d\n", g_data.match_count);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}
