#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <string.h>

// --- Mersenne Twister (Do Not Modify - Include This Verbatim) ---
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

// --- Benchmark Globals ---
int sequence_length;
char *rna_sequence;
int **dp_table;
int final_result;

// --- Helper Function ---
// Determines if two RNA bases can form a Watson-Crick pair.
int can_pair(char c1, char c2) {
    return ((c1 == 'A' && c2 == 'U') || (c1 == 'U' && c2 == 'A') ||
            (c1 == 'G' && c2 == 'C') || (c1 == 'C' && c2 == 'G'));
}

// --- Benchmark Functions ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <sequence_length> <seed>\n", argv[0]);
        exit(1);
    }

    sequence_length = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if (sequence_length <= 0) {
        fprintf(stderr, "FATAL: sequence_length must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate memory for the RNA sequence
    rna_sequence = (char *)malloc(sequence_length * sizeof(char));
    if (!rna_sequence) {
        fprintf(stderr, "FATAL: Failed to allocate memory for RNA sequence.\n");
        exit(1);
    }
    
    // Generate a random RNA sequence (A, C, G, U)
    const char bases[] = {'A', 'C', 'G', 'U'};
    for (int i = 0; i < sequence_length; ++i) {
        rna_sequence[i] = bases[mt_rand() % 4];
    }
    
    // Allocate the dynamic programming table
    dp_table = (int **)malloc(sequence_length * sizeof(int *));
    if (!dp_table) {
        fprintf(stderr, "FATAL: Failed to allocate memory for DP table rows.\n");
        free(rna_sequence);
        exit(1);
    }
    for (int i = 0; i < sequence_length; ++i) {
        dp_table[i] = (int *)malloc(sequence_length * sizeof(int));
        if (!dp_table[i]) {
            fprintf(stderr, "FATAL: Failed to allocate memory for DP table column %d.\n", i);
            // Cleanup already allocated memory
            for (int k = 0; k < i; ++k) free(dp_table[k]);
            free(dp_table);
            free(rna_sequence);
            exit(1);
        }
        // Initialize score to 0. Required by the DP algorithm.
        memset(dp_table[i], 0, sequence_length * sizeof(int));
    }
}

void run_computation() {
    // Implementation of the Nussinov algorithm for RNA secondary structure prediction.
    // This algorithm maximizes the number of base pairs in a given RNA sequence.
    // The complexity is O(N^3) due to the three nested loops (len, i, k).
    // dp_table[i][j] stores the maximum number of pairs in the subsequence from i to j.

    for (int len = 1; len < sequence_length; ++len) {
        for (int i = 0; i < sequence_length - len; ++i) {
            int j = i + len;

            // Case 1: Base j is unpaired. The score is taken from the shorter subsequence.
            int score_unpaired = (j > 0) ? dp_table[i][j - 1] : 0;

            // Case 2: Base j pairs with some base k where i <= k < j.
            int max_bifurcation_score = 0;
            for (int k = i; k < j; ++k) {
                if (can_pair(rna_sequence[k], rna_sequence[j])) {
                    // Score contributed by the segment before the pair (if any)
                    int score_before = (k > i) ? dp_table[i][k - 1] : 0;
                    // Score contributed by the segment enclosed by the pair k,j (if any)
                    int score_enclosed = (k + 1 < j) ? dp_table[k + 1][j - 1] : 0;
                    
                    int current_bifurcation_score = score_before + score_enclosed + 1; // +1 for the new pair (k,j)
                    if (current_bifurcation_score > max_bifurcation_score) {
                        max_bifurcation_score = current_bifurcation_score;
                    }
                }
            }
            
            // The final score for dp[i][j] is the max of j being unpaired vs. j being paired.
            dp_table[i][j] = (score_unpaired > max_bifurcation_score) ? score_unpaired : max_bifurcation_score;
        }
    }
    
    // The final result is the max pairs for the entire sequence.
    final_result = dp_table[0][sequence_length - 1];
}

void cleanup() {
    for (int i = 0; i < sequence_length; ++i) {
        free(dp_table[i]);
    }
    free(dp_table);
    free(rna_sequence);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final result (max number of pairs) to stdout
    printf("%d\n", final_result);

    // Print the execution time to stderr
    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}
