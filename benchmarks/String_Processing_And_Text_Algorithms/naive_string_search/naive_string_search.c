#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// Mersenne Twister (DO NOT MODIFY)
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
// End of Mersenne Twister

// --- Benchmark Globals ---
struct {
    char *text;
    char *pattern;
    size_t text_len;
    size_t pattern_len;
    int total_matches;
} g_data;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <text_length> <pattern_length> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.text_len = strtoul(argv[1], NULL, 10);
    g_data.pattern_len = strtoul(argv[2], NULL, 10);
    uint32_t seed = strtoul(argv[3], NULL, 10);
    
    if (g_data.text_len == 0 || g_data.pattern_len == 0) {
        fprintf(stderr, "Error: text and pattern lengths must be greater than 0.\n");
        exit(1);
    }
    if (g_data.pattern_len > g_data.text_len) {
        fprintf(stderr, "Error: pattern length cannot be greater than text length.\n");
        exit(1);
    }

    mt_seed(seed);

    g_data.text = (char *)malloc(g_data.text_len + 1);
    g_data.pattern = (char *)malloc(g_data.pattern_len + 1);

    if (!g_data.text || !g_data.pattern) {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        exit(1);
    }

    // Generate random printable ASCII characters (' ' to '~')
    for (size_t i = 0; i < g_data.text_len; i++) {
        g_data.text[i] = (char)(mt_rand() % 95 + 32);
    }
    g_data.text[g_data.text_len] = '\0';

    for (size_t i = 0; i < g_data.pattern_len; i++) {
        g_data.pattern[i] = (char)(mt_rand() % 95 + 32);
    }
    g_data.pattern[g_data.pattern_len] = '\0';

    g_data.total_matches = 0;
}

void run_computation() {
    int count = 0;
    size_t n = g_data.text_len;
    size_t m = g_data.pattern_len;

    // Naive string search algorithm
    for (size_t i = 0; i <= n - m; i++) {
        int match = 1;
        for (size_t j = 0; j < m; j++) {
            if (g_data.text[i + j] != g_data.pattern[j]) {
                match = 0;
                break;
            }
        }
        if (match) {
            count++;
        }
    }
    g_data.total_matches = count;
}

void cleanup() {
    free(g_data.text);
    free(g_data.pattern);
}

// --- Main Execution Logic ---

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
    // Print result to stdout
    printf("%d\n", g_data.total_matches);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
