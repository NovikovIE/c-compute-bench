#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- Mersenne Twister (MT19937) random number generator ---
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

// --- Benchmark-specific code ---

// Parameters and data structures
static size_t SEQ1_LEN;
static size_t SEQ2_LEN;
static char *seq1;
static char *seq2;
static int *score_matrix;

// Scoring constants
#define MATCH_SCORE 5
#define MISMATCH_SCORE -3
#define GAP_PENALTY -4

// Result of the computation
static int final_max_score;

#define MAX(a, b) ((a) > (b) ? (a) : (b))

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <sequence1_length> <sequence2_length> <seed>\n", argv[0]);
        exit(1);
    }

    SEQ1_LEN = (size_t)atoi(argv[1]);
    SEQ2_LEN = (size_t)atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (SEQ1_LEN == 0 || SEQ2_LEN == 0) {
        fprintf(stderr, "FATAL: Sequence lengths must be positive integers.\n");
        exit(1);
    }

    mt_seed(seed);

    const char alphabet[] = "ACGT";

    seq1 = (char *)malloc(SEQ1_LEN * sizeof(char));
    seq2 = (char *)malloc(SEQ2_LEN * sizeof(char));
    // Use calloc to initialize the scoring matrix to all zeros, as required by the algorithm.
    score_matrix = (int *)calloc((SEQ1_LEN + 1) * (SEQ2_LEN + 1), sizeof(int));

    if (!seq1 || !seq2 || !score_matrix) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Generate random DNA sequences
    for (size_t i = 0; i < SEQ1_LEN; ++i) {
        seq1[i] = alphabet[mt_rand() % 4];
    }
    for (size_t i = 0; i < SEQ2_LEN; ++i) {
        seq2[i] = alphabet[mt_rand() % 4];
    }
}

void run_computation() {
    int max_score_so_far = 0;
    
    for (size_t i = 1; i <= SEQ1_LEN; ++i) {
        for (size_t j = 1; j <= SEQ2_LEN; ++j) {
            int score_sub = (seq1[i - 1] == seq2[j - 1]) ? MATCH_SCORE : MISMATCH_SCORE;
            
            int diag_score = score_matrix[(i - 1) * (SEQ2_LEN + 1) + (j - 1)] + score_sub;
            int del_score = score_matrix[(i - 1) * (SEQ2_LEN + 1) + j] + GAP_PENALTY;
            int ins_score = score_matrix[i * (SEQ2_LEN + 1) + (j - 1)] + GAP_PENALTY;
            
            int current_score = MAX(0, MAX(diag_score, MAX(del_score, ins_score)));
            
            score_matrix[i * (SEQ2_LEN + 1) + j] = current_score;

            if (current_score > max_score_so_far) {
                max_score_so_far = current_score;
            }
        }
    }
    final_max_score = max_score_so_far;
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

    // Print the final accumulated result to stdout
    printf("%d\n", final_max_score);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
