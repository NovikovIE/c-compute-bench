#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <string.h>

// --- Mersenne Twister (MT19937) Generator (DO NOT MODIFY) ---
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

// Benchmark parameters and data
size_t sequence1_length;
size_t sequence2_length;
char *seq1;
char *seq2;
int *score_matrix;

// Result
int final_score;

// Scoring constants for the alignment algorithm
#define MATCH_SCORE 1
#define MISMATCH_SCORE -1
#define GAP_PENALTY -2

// Helper function for finding the maximum of three integers
static inline int max3(int a, int b, int c) {
    int max_ab = (a > b) ? a : b;
    return (max_ab > c) ? max_ab : c;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <sequence1_length> <sequence2_length> <seed>\n", argv[0]);
        exit(1);
    }

    sequence1_length = atol(argv[1]);
    sequence2_length = atol(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);
    mt_seed(seed);

    // Allocate memory for the two sequences
    // Add 1 for a potential null terminator, although it's not strictly required by the algorithm itself
    seq1 = (char *)malloc((sequence1_length + 1) * sizeof(char));
    seq2 = (char *)malloc((sequence2_length + 1) * sizeof(char));
    if (!seq1 || !seq2) {
        fprintf(stderr, "FATAL: Memory allocation failed for sequences.\n");
        exit(1);
    }

    // Generate random DNA sequences using a simple 4-character alphabet (A, C, G, T)
    const char alphabet[] = "ACGT";
    for (size_t i = 0; i < sequence1_length; ++i) {
        seq1[i] = alphabet[mt_rand() % 4];
    }
    seq1[sequence1_length] = '\0';
    for (size_t i = 0; i < sequence2_length; ++i) {
        seq2[i] = alphabet[mt_rand() % 4];
    }
    seq2[sequence2_length] = '\0';

    // Allocate memory for the score matrix. 
    // The dimensions are (sequence1_length + 1) x (sequence2_length + 1)
    size_t matrix_rows = sequence1_length + 1;
    size_t matrix_cols = sequence2_length + 1;
    score_matrix = (int *)malloc(matrix_rows * matrix_cols * sizeof(int));
    if (!score_matrix) {
        fprintf(stderr, "FATAL: Memory allocation for score matrix failed.\n");
        free(seq1);
        free(seq2);
        exit(1);
    }
}

void run_computation() {
    size_t m = sequence1_length;
    size_t n = sequence2_length;
    size_t cols = n + 1;

    // --- Initialization Step ---
    // Initialize the first row and column of the matrix with gap penalties
    score_matrix[0] = 0;
    for (size_t i = 1; i <= m; ++i) {
        score_matrix[i * cols] = i * GAP_PENALTY;
    }
    for (size_t j = 1; j <= n; ++j) {
        score_matrix[j] = j * GAP_PENALTY;
    }

    // --- Matrix Filling Step ---
    // Iterate through the matrix, filling each cell with the optimal score
    for (size_t i = 1; i <= m; ++i) {
        for (size_t j = 1; j <= n; ++j) {
            // Score if characters align (match or mismatch)
            int match_val = (seq1[i - 1] == seq2[j - 1]) ? MATCH_SCORE : MISMATCH_SCORE;
            int score_diag = score_matrix[(i - 1) * cols + (j - 1)] + match_val;

            // Score if introducing a gap in sequence2 (moving down)
            int score_up = score_matrix[(i - 1) * cols + j] + GAP_PENALTY;

            // Score if introducing a gap in sequence1 (moving right)
            int score_left = score_matrix[i * cols + (j - 1)] + GAP_PENALTY;

            // The cell's score is the maximum of the three possibilities
            score_matrix[i * cols + j] = max3(score_diag, score_up, score_left);
        }
    }

    // The final alignment score is in the bottom-right corner of the matrix
    final_score = score_matrix[m * cols + n];
}

void cleanup() {
    free(seq1);
    free(seq2);
    free(score_matrix);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final alignment score to stdout to prevent dead code elimination
    printf("%d\n", final_score);

    // Print the measured time to stderr for parsing
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
