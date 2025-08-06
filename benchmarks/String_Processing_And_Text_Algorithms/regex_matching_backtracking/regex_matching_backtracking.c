#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>

// --- BEGIN MERSENNE TWISTER (DO NOT MODIFY) ---
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
// --- END MERSENNE TWISTER ---

// --- BENCHMARK DATA ---
char* text;
char* regex_pattern;
long long total_states_explored; // Use as result to prevent dead code elimination

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <text_length> <regex_complexity_score> <seed>\n", argv[0]);
        exit(1);
    }

    int text_length = atoi(argv[1]);
    int regex_complexity_score = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    mt_seed(seed);

    total_states_explored = 0;

    // Allocate memory for the text
    text = (char*)malloc((text_length + 1) * sizeof(char));
    if (!text) {
        fprintf(stderr, "Failed to allocate memory for text.\n");
        exit(1);
    }

    // Generate a simple text that is a near-miss for the regex pattern
    // This forces the backtracking engine to explore all paths before failing.
    for (int i = 0; i < text_length; i++) {
        text[i] = 'a';
    }
    text[text_length] = '\0';

    // Allocate memory for the regex pattern
    // Pattern is of the form: a*a*a*...b (complexity_score determines number of 'a*' pairs)
    int pattern_len = regex_complexity_score * 2 + 1; // e.g., score 3 -> a*a*a*b
    regex_pattern = (char*)malloc((pattern_len + 1) * sizeof(char));
    if (!regex_pattern) {
        fprintf(stderr, "Failed to allocate memory for regex pattern.\n");
        exit(1);
    }

    // Generate the regex pattern that causes catastrophic backtracking
    char* p = regex_pattern;
    for (int i = 0; i < regex_complexity_score; i++) {
        *p++ = 'a';
        *p++ = '*';
    }
    *p++ = 'b'; // The final character that will never match
    *p = '\0';
}

// A simple recursive backtracking regex matcher supporting '.', and '*'.
// This function's structure is designed to exhibit exponential runtime on specific patterns.
static int match_recursive(const char* regex, const char* text) {
    total_states_explored++;

    if (regex[0] == '\0') {
        return 1; // End of regex, successful match.
    }

    if (regex[1] == '*') {
        // Handle the '*' quantifier: non-deterministic choice
        // Path 1: Match zero characters and skip the 'c*' pattern.
        if (match_recursive(regex + 2, text)) {
            return 1;
        }
        // Path 2: If Path 1 fails, match one character and repeat the 'c*' pattern.
        if (text[0] != '\0' && (regex[0] == '.' || regex[0] == text[0])) {
            return match_recursive(regex, text + 1);
        }
    }

    // Standard character match
    if (text[0] != '\0' && (regex[0] == '.' || regex[0] == text[0])) {
        return match_recursive(regex + 1, text + 1);
    }

    return 0; // No match found
}

void run_computation() {
    // The match will always fail, but the computation consists of exploring
    // the entire state space, which is counted in 'total_states_explored'.
    match_recursive(regex_pattern, text);
}

void cleanup() {
    free(text);
    free(regex_pattern);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the accumulated result to stdout to prevent dead code elimination.
    printf("%lld\n", total_states_explored);

    // Print the timing information to stderr.
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
