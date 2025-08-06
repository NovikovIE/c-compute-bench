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

// --- Benchmark Globals ---
static char* text;
static char* pattern;
static size_t text_length;
static size_t pattern_length;
static long long final_result = 0;

// --- Benchmark Functions ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <text_length> <pattern_length> <seed>\n", argv[0]);
        exit(1);
    }

    text_length = strtoul(argv[1], NULL, 10);
    pattern_length = strtoul(argv[2], NULL, 10);
    uint32_t seed = (uint32_t)strtoul(argv[3], NULL, 10);

    if (text_length == 0 || pattern_length == 0) {
        fprintf(stderr, "FATAL: text and pattern lengths must be positive.\n");
        exit(1);
    }
    if (pattern_length > text_length) {
        fprintf(stderr, "FATAL: pattern length cannot be greater than text length.\n");
        exit(1);
    }

    mt_seed(seed);

    text = (char*)malloc(text_length * sizeof(char));
    if (!text) {
        fprintf(stderr, "FATAL: Memory allocation failed for text.\n");
        exit(1);
    }

    pattern = (char*)malloc(pattern_length * sizeof(char));
    if (!pattern) {
        free(text);
        fprintf(stderr, "FATAL: Memory allocation failed for pattern.\n");
        exit(1);
    }

    // Generate random character strings (e.g., lowercase alphabet)
    for (size_t i = 0; i < text_length; ++i) {
        text[i] = 'a' + (mt_rand() % 26);
    }

    for (size_t i = 0; i < pattern_length; ++i) {
        pattern[i] = 'a' + (mt_rand() % 26);
    }

    final_result = 0;
}

void cleanup() {
    free(text);
    free(pattern);
}

// Helper to compute Z-array for a string str.
static void compute_z_array(const char *str, int *z, size_t n) {
    size_t l = 0, r = 0;
    z[0] = 0; // Z[0] is not used in pattern matching

    for (size_t i = 1; i < n; ++i) {
        if (i > r) {
            // Case 1: i is outside the current [l, r] Z-box.
            // We need to compute Z[i] naively.
            l = r = i;
            while (r < n && str[r - l] == str[r]) {
                r++;
            }
            z[i] = r - l;
            if (r > 0) r--; // r is the inclusive right boundary
        } else {
            // Case 2: i is inside the current Z-box.
            size_t k = i - l;
            if (z[k] < r - i + 1) {
                // Z[k] value does not extend to or beyond r.
                z[i] = z[k];
            } else {
                // Z[k] value extends to or beyond r.
                // We need to start a new match from r.
                l = i;
                while (r < n && str[r - l] == str[r]) {
                    r++;
                }
                z[i] = r - l;
                if(r > 0) r--;
            }
        }
    }
}

void run_computation() {
    size_t concat_len = pattern_length + 1 + text_length;

    char *concat_str = (char *)malloc(concat_len * sizeof(char));
    if (!concat_str) return; // Fail gracefully

    int *z_array = (int *)malloc(concat_len * sizeof(int));
    if (!z_array) {
        free(concat_str);
        return; // Fail gracefully
    }

    // Create the concatenated string: pattern + '$' + text
    memcpy(concat_str, pattern, pattern_length);
    concat_str[pattern_length] = '$'; // A separator not in the alphabet
    memcpy(concat_str + pattern_length + 1, text, text_length);

    // Compute the Z-array for the concatenated string
    compute_z_array(concat_str, z_array, concat_len);

    // Find matches by checking the Z-array
    long long sum_indices = 0;
    for (size_t i = pattern_length + 1; i < concat_len; ++i) {
        if (z_array[i] == pattern_length) {
            // Match found at index (i - pattern_length - 1) in the original text
            sum_indices += (i - pattern_length - 1);
        }
    }
    final_result = sum_indices;

    free(concat_str);
    free(z_array);
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
    printf("%lld\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
