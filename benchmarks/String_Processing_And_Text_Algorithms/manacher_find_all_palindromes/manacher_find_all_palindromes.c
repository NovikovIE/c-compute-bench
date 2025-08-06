#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

#define MIN(a, b) (((a) < (b)) ? (a) : (b))

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
// --- End of MT19937 ---

// --- Benchmark Data and Globals ---
typedef struct {
    int text_length;
    char* text;
    char* transformed_text; // For Manacher's algorithm
    int* p_array;           // Stores palindrome lengths
    long long total_palindromes;
} BenchmarkData;

BenchmarkData g_data;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <text_length> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.text_length = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);
    mt_seed(seed);

    if (g_data.text_length <= 0) {
        fprintf(stderr, "Error: text_length must be a positive integer.\n");
        exit(1);
    }

    // 1. Allocate and generate the original text
    g_data.text = (char*)malloc((g_data.text_length + 1) * sizeof(char));
    if (!g_data.text) {
        fprintf(stderr, "Memory allocation failed for text.\n");
        exit(1);
    }
    const char alphabet[] = "abcdefghijklmnopqrstuvwxyz";
    int alphabet_size = sizeof(alphabet) - 1;
    for (int i = 0; i < g_data.text_length; ++i) {
        g_data.text[i] = alphabet[mt_rand() % alphabet_size];
    }
    g_data.text[g_data.text_length] = '\0';

    // 2. Allocate and create the transformed text for Manacher's algorithm
    // Transformed text is of form: ^#c1#c2#...#cn#$
    int transformed_len = 2 * g_data.text_length + 3;
    g_data.transformed_text = (char*)malloc((transformed_len + 1) * sizeof(char));
    if (!g_data.transformed_text) {
        fprintf(stderr, "Memory allocation failed for transformed_text.\n");
        exit(1);
    }
    
    g_data.transformed_text[0] = '^'; // Start sentinel
    for (int i = 0; i < g_data.text_length; ++i) {
        g_data.transformed_text[2 * i + 1] = '#';
        g_data.transformed_text[2 * i + 2] = g_data.text[i];
    }
    g_data.transformed_text[transformed_len - 2] = '#';
    g_data.transformed_text[transformed_len - 1] = '$'; // End sentinel
    g_data.transformed_text[transformed_len] = '\0';

    // 3. Allocate the palindrome length array (P array)
    g_data.p_array = (int*)malloc(transformed_len * sizeof(int));
    if (!g_data.p_array) {
        fprintf(stderr, "Memory allocation failed for p_array.\n");
        exit(1);
    }
    // No need to initialize p_array, it will be filled by the algorithm

    // 4. Initialize result accumulator
    g_data.total_palindromes = 0;
}

void run_computation() {
    int n = g_data.text_length;
    int transformed_len = 2 * n + 3;
    char* S = g_data.transformed_text;
    int* P = g_data.p_array;

    int center = 0, right_boundary = 0;

    for (int i = 1; i < transformed_len - 1; i++) {
        int i_mirror = 2 * center - i;

        // If i is within the current right_boundary, we can leverage the palindrome length
        // from its mirror position. Otherwise, start with a length of 0.
        P[i] = (right_boundary > i) ? MIN(right_boundary - i, P[i_mirror]) : 0;

        // Attempt to expand palindrome centered at i.
        // The sentinels '^' and '$' prevent running off the ends of the string.
        while (S[i + 1 + P[i]] == S[i - 1 - P[i]]) {
            P[i]++;
        }

        // If palindrome centered at i expands past right_boundary, update center and right_boundary.
        if (i + P[i] > right_boundary) {
            center = i;
            right_boundary = i + P[i];
        }
    }

    // Sum up the number of palindromes.
    // The length of the palindrome in the original string corresponding to P[i] is P[i].
    // The number of palindromic substrings centered at this point is (P[i] + 1) / 2.
    long long count = 0;
    for (int i = 1; i < transformed_len - 1; i++) {
        count += (long long)(P[i] + 1) / 2;
    }
    g_data.total_palindromes = count;
}

void cleanup() {
    free(g_data.text);
    free(g_data.transformed_text);
    free(g_data.p_array);
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
    printf("%lld\n", g_data.total_palindromes);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
