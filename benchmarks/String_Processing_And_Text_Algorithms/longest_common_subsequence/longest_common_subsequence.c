#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

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
// --- End of Mersenne Twister ---

// --- Benchmark Globals ---
static int m, n;
static char* string1;
static char* string2;
static int* dp_table;
static int lcs_length_result;

// --- Utility Macro ---
#define MAX(a, b) ((a) > (b) ? (a) : (b))

// --- Benchmark Functions ---
void setup_benchmark(int argc, char* argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <string1_length> <string2_length> <seed>\n", argv[0]);
        exit(1);
    }

    m = atoi(argv[1]);
    n = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (m <= 0 || n <= 0) {
        fprintf(stderr, "ERROR: String lengths must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    string1 = (char*)malloc((size_t)m * sizeof(char));
    string2 = (char*)malloc((size_t)n * sizeof(char));
    
    size_t dp_table_size = (size_t)(m + 1) * (n + 1);
    dp_table = (int*)calloc(dp_table_size, sizeof(int));

    if (!string1 || !string2 || !dp_table) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        free(string1);
        free(string2);
        free(dp_table);
        exit(1);
    }

    // Using a 4-character alphabet for simplicity (e.g., DNA)
    const char alphabet[] = "ACGT";
    for (int i = 0; i < m; ++i) {
        string1[i] = alphabet[mt_rand() % 4];
    }
    for (int i = 0; i < n; ++i) {
        string2[i] = alphabet[mt_rand() % 4];
    }
}

void run_computation() {
    size_t width = n + 1;
    // The DP table is already initialized to 0s by calloc,
    // which covers the base cases where i=0 or j=0.
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            if (string1[i - 1] == string2[j - 1]) {
                dp_table[i * width + j] = 1 + dp_table[(i - 1) * width + (j - 1)];
            } else {
                dp_table[i * width + j] = MAX(dp_table[(i - 1) * width + j], dp_table[i * width + (j - 1)]);
            }
        }
    }
    lcs_length_result = dp_table[m * width + n];
}

void cleanup() {
    free(string1);
    free(string2);
    free(dp_table);
}

// --- Main Function ---
int main(int argc, char* argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final result to stdout
    printf("%d\n", lcs_length_result);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
