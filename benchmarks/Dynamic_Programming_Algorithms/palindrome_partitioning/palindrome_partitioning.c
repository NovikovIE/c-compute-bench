/*
 * BENCHMARK: Palindrome Partitioning
 * DESCRIPTION: This benchmark calculates the minimum number of cuts needed to partition a 
 *              string into palindromic substrings. It uses a dynamic programming approach.
 *              First, it precomputes a table identifying all palindromic substrings.
 *              Then, it uses this table to determine the minimum cuts for prefixes of
 *              increasing length.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

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

// --- BENCHMARK DATA AND FUNCTIONS ---

typedef struct {
    int string_length;
    char *input_string;
    int **is_palindrome; // Using int for boolean (0/1)
    int *cuts;
    int final_result;
} BenchmarkData;

BenchmarkData g_data;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <string_length> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.string_length = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if (g_data.string_length <= 0) {
        fprintf(stderr, "FATAL: string_length must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    g_data.input_string = (char *)malloc((g_data.string_length + 1) * sizeof(char));
    if (!g_data.input_string) {
        fprintf(stderr, "FATAL: Memory allocation failed for input_string.\n");
        exit(1);
    }

    for (int i = 0; i < g_data.string_length; ++i) {
        g_data.input_string[i] = 'a' + (mt_rand() % 26);
    }
    g_data.input_string[g_data.string_length] = '\0';

    g_data.is_palindrome = (int **)malloc(g_data.string_length * sizeof(int *));
    if (!g_data.is_palindrome) {
        fprintf(stderr, "FATAL: Memory allocation failed for is_palindrome pointers.\n");
        exit(1);
    }
    for (int i = 0; i < g_data.string_length; ++i) {
        g_data.is_palindrome[i] = (int *)malloc(g_data.string_length * sizeof(int));
        if (!g_data.is_palindrome[i]) {
            fprintf(stderr, "FATAL: Memory allocation failed for is_palindrome row %d.\n", i);
            exit(1);
        }
    }

    g_data.cuts = (int *)malloc(g_data.string_length * sizeof(int));
    if (!g_data.cuts) {
        fprintf(stderr, "FATAL: Memory allocation failed for cuts array.\n");
        exit(1);
    }

    g_data.final_result = 0;
}

void run_computation() {
    int n = g_data.string_length;
    if (n == 0) {
        g_data.final_result = 0;
        return;
    }

    // Step 1: Precompute the is_palindrome table.
    // is_palindrome[i][j] is true if the substring s[i..j] is a palindrome.
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            g_data.is_palindrome[i][j] = 0;
        }
    }

    for (int len = 1; len <= n; ++len) {
        for (int i = 0; i <= n - len; ++i) {
            int j = i + len - 1;
            if (len == 1) {
                g_data.is_palindrome[i][j] = 1;
            } else if (len == 2) {
                g_data.is_palindrome[i][j] = (g_data.input_string[i] == g_data.input_string[j]);
            } else {
                g_data.is_palindrome[i][j] = (g_data.input_string[i] == g_data.input_string[j]) && g_data.is_palindrome[i + 1][j - 1];
            }
        }
    }

    // Step 2: Compute the minimum cuts.
    // cuts[i] is the minimum number of cuts for the prefix s[0..i].
    for (int i = 0; i < n; ++i) {
        if (g_data.is_palindrome[0][i]) {
            g_data.cuts[i] = 0;
        } else {
            g_data.cuts[i] = i; // Max possible cuts
            for (int j = 1; j <= i; ++j) {
                if (g_data.is_palindrome[j][i]) {
                    int new_cuts = g_data.cuts[j - 1] + 1;
                    if (new_cuts < g_data.cuts[i]) {
                        g_data.cuts[i] = new_cuts;
                    }
                }
            }
        }
    }

    g_data.final_result = g_data.cuts[n - 1];
}

void cleanup() {
    if (g_data.is_palindrome) {
        for (int i = 0; i < g_data.string_length; ++i) {
            free(g_data.is_palindrome[i]);
        }
        free(g_data.is_palindrome);
    }
    free(g_data.input_string);
    free(g_data.cuts);
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
    printf("%d\n", g_data.final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
