#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

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

// --- Benchmark Data Structures ---
typedef struct {
    char *pattern;
    int token_id;
    int pattern_length;
} LexicalRule;

typedef struct {
    size_t input_size;
    int num_rules;
    char* input_string;
    LexicalRule* rules;
    long long total_tokens_value;
} BenchmarkData;

BenchmarkData g_data;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <input_size_bytes> <num_lexical_rules> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.input_size = atol(argv[1]);
    g_data.num_rules = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    mt_seed(seed);

    // Allocate and generate input string
    g_data.input_string = (char*)malloc(g_data.input_size * sizeof(char));
    if (!g_data.input_string) {
        fprintf(stderr, "Failed to allocate memory for input string.\n");
        exit(1);
    }
    for (size_t i = 0; i < g_data.input_size; ++i) {
        // Use a small alphabet to increase match probability
        g_data.input_string[i] = 'a' + (mt_rand() % 8);
    }

    // Allocate and generate lexical rules
    g_data.rules = (LexicalRule*)malloc(g_data.num_rules * sizeof(LexicalRule));
    if (!g_data.rules) {
        fprintf(stderr, "Failed to allocate memory for rules.\n");
        exit(1);
    }

    for (int i = 0; i < g_data.num_rules; ++i) {
        int len = 2 + (mt_rand() % 7); // Pattern length between 2 and 8
        g_data.rules[i].pattern_length = len;
        g_data.rules[i].token_id = mt_rand() % 1000;
        g_data.rules[i].pattern = (char*)malloc((len + 1) * sizeof(char));
        if (!g_data.rules[i].pattern) {
            fprintf(stderr, "Failed to allocate memory for a rule pattern.\n");
            exit(1);
        }

        for (int j = 0; j < len; ++j) {
            g_data.rules[i].pattern[j] = 'a' + (mt_rand() % 8);
        }
        g_data.rules[i].pattern[len] = '\0';
    }
    g_data.total_tokens_value = 0;
}

void run_computation() {
    long long accumulated_value = 0;
    size_t current_pos = 0;

    while (current_pos < g_data.input_size) {
        int best_match_len = 0;
        int best_match_token_id = -1;

        // Find the longest matching rule at the current position
        for (int i = 0; i < g_data.num_rules; ++i) {
            int p_len = g_data.rules[i].pattern_length;
            if (p_len > best_match_len && (current_pos + p_len) <= g_data.input_size) {
                if (strncmp(g_data.input_string + current_pos, g_data.rules[i].pattern, p_len) == 0) {
                    best_match_len = p_len;
                    best_match_token_id = g_data.rules[i].token_id;
                }
            }
        }

        if (best_match_len > 0) {
            accumulated_value += best_match_token_id;
            current_pos += best_match_len;
        } else {
            current_pos++;
        }
    }
    g_data.total_tokens_value = accumulated_value;
}

void cleanup() {
    if (g_data.rules) {
        for (int i = 0; i < g_data.num_rules; i++) {
            free(g_data.rules[i].pattern);
        }
        free(g_data.rules);
    }
    free(g_data.input_string);
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
    printf("%lld\n", g_data.total_tokens_value);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
