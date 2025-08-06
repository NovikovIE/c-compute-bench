#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- Mersenne Twister (DO NOT MODIFY) ---
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
// --- End of Mersenne Twister ---

// --- Benchmark Globals ---
static int string1_length;
static int string2_length;
static char *s1;
static char *s2;
static int *dp_matrix; // DP matrix, stored as a 1D array
static int final_result;

// --- Helper Functions ---
#define MIN3(a, b, c) ((a) < (b) ? ((a) < (c) ? (a) : (c)) : ((b) < (c) ? (b) : (c)))

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <string1_length> <string2_length> <seed>\n", argv[0]);
        exit(1);
    }

    string1_length = atoi(argv[1]);
    string2_length = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (string1_length <= 0 || string2_length <= 0) {
        fprintf(stderr, "FATAL: String lengths must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate memory for the two strings (+1 for null terminator)
    s1 = (char *)malloc((string1_length + 1) * sizeof(char));
    s2 = (char *)malloc((string2_length + 1) * sizeof(char));
    if (!s1 || !s2) {
        fprintf(stderr, "FATAL: Memory allocation failed for strings.\n");
        exit(1);
    }

    // Generate random strings from a small alphabet
    for (int i = 0; i < string1_length; i++) {
        s1[i] = 'a' + (mt_rand() % 26);
    }
    s1[string1_length] = '\0';

    for (int i = 0; i < string2_length; i++) {
        s2[i] = 'a' + (mt_rand() % 26);
    }
    s2[string2_length] = '\0';

    // Allocate memory for the DP matrix
    // Size is (string1_length+1) x (string2_length+1)
    size_t dp_size = (size_t)(string1_length + 1) * (string2_length + 1);
    dp_matrix = (int *)malloc(dp_size * sizeof(int));
    if (!dp_matrix) {
        fprintf(stderr, "FATAL: Memory allocation failed for DP matrix.\n");
        free(s1);
        free(s2);
        exit(1);
    }
}

void run_computation() {
    int n2_plus_1 = string2_length + 1;

    // Initialize the first row and column of the DP matrix
    for (int i = 0; i <= string1_length; i++) {
        dp_matrix[i * n2_plus_1 + 0] = i;
    }
    for (int j = 0; j <= string2_length; j++) {
        dp_matrix[0 * n2_plus_1 + j] = j;
    }

    // Fill the rest of the DP matrix
    for (int i = 1; i <= string1_length; i++) {
        for (int j = 1; j <= string2_length; j++) {
            int cost = (s1[i - 1] == s2[j - 1]) ? 0 : 1;
            int deletion_cost = dp_matrix[(i - 1) * n2_plus_1 + j] + 1;
            int insertion_cost = dp_matrix[i * n2_plus_1 + (j - 1)] + 1;
            int substitution_cost = dp_matrix[(i - 1) * n2_plus_1 + (j - 1)] + cost;

            dp_matrix[i * n2_plus_1 + j] = MIN3(deletion_cost, insertion_cost, substitution_cost);
        }
    }

    final_result = dp_matrix[string1_length * n2_plus_1 + string2_length];
}

void cleanup() {
    free(s1);
    free(s2);
    free(dp_matrix);
}

// --- Main Function ---

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%d\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
