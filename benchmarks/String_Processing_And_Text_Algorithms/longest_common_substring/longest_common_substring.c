#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- BEGIN MERSENNE TWISTER --- 
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

// Global variables for data sharing
char *string1 = NULL;
char *string2 = NULL;
int **dp_table = NULL;
size_t S1_LEN, S2_LEN;
int final_result_length = 0;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <string1_length> <string2_length> <seed>\n", argv[0]);
        exit(1);
    }

    S1_LEN = atoi(argv[1]);
    S2_LEN = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    mt_seed(seed);

    // Allocate memory for strings (+1 for null terminator)
    string1 = (char *)malloc((S1_LEN + 1) * sizeof(char));
    string2 = (char *)malloc((S2_LEN + 1) * sizeof(char));
    if (!string1 || !string2) {
        fprintf(stderr, "Failed to allocate memory for strings.\n");
        exit(1);
    }

    // Generate random strings from a small alphabet to increase common substrings
    for (size_t i = 0; i < S1_LEN; ++i) {
        string1[i] = 'A' + (mt_rand() % 4); // Alphabet: A, B, C, D
    }
    string1[S1_LEN] = '\0';

    for (size_t i = 0; i < S2_LEN; ++i) {
        string2[i] = 'A' + (mt_rand() % 4);
    }
    string2[S2_LEN] = '\0';

    // Allocate DP table for computation phase, initialized to zero
    dp_table = (int **)malloc((S1_LEN + 1) * sizeof(int *));
    if (!dp_table) {
        fprintf(stderr, "Failed to allocate memory for DP table rows.\n");
        exit(1);
    }
    for (size_t i = 0; i <= S1_LEN; ++i) {
        dp_table[i] = (int *)calloc(S2_LEN + 1, sizeof(int));
        if (!dp_table[i]) {
             fprintf(stderr, "Failed to allocate memory for DP table columns.\n");
             exit(1);
        }
    }
}

void run_computation() {
    int max_len = 0;
    // Dynamic programming to find the length of the longest common substring
    // dp_table[i][j] contains length of longest common suffix of string1[0..i-1] and string2[0..j-1]
    for (size_t i = 1; i <= S1_LEN; i++) {
        for (size_t j = 1; j <= S2_LEN; j++) {
            if (string1[i - 1] == string2[j - 1]) {
                dp_table[i][j] = dp_table[i - 1][j - 1] + 1;
                if (dp_table[i][j] > max_len) {
                    max_len = dp_table[i][j];
                }
            } else {
                dp_table[i][j] = 0;
            }
        }
    }
    final_result_length = max_len;
}

void cleanup() {
    if (dp_table) {
        for (size_t i = 0; i <= S1_LEN; ++i) {
            free(dp_table[i]);
        }
        free(dp_table);
    }
    free(string1);
    free(string2);
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
    printf("%d\n", final_result_length);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
