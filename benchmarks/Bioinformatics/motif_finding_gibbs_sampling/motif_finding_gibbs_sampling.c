#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- BEGIN MERSENNE TWISTER --- (Do Not Modify)
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
// --- END MERSENNE TWISTER ---

// Global struct for benchmark data and parameters
typedef struct {
    // Parameters
    int num_sequences;
    int average_sequence_length;
    int motif_length;
    int num_iterations;
    uint32_t seed;

    // Data structures
    char** sequences;
    int* sequence_lengths;
    int* motif_starts; 
    double** profile_matrix; // Position Weight Matrix (PWM)
    double* score_buffer; 

    // Final result
    long long final_result;
} BenchmarkData;

BenchmarkData g_data;

static inline int nucleotide_to_index(char n) {
    switch (n) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return 0; // Should not happen
    }
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s num_sequences avg_seq_len motif_len num_iter seed\n", argv[0]);
        exit(1);
    }

    g_data.num_sequences = atoi(argv[1]);
    g_data.average_sequence_length = atoi(argv[2]);
    g_data.motif_length = atoi(argv[3]);
    g_data.num_iterations = atoi(argv[4]);
    g_data.seed = (uint32_t)atoi(argv[5]);

    mt_seed(g_data.seed);

    // Allocate memory for sequences and their lengths
    g_data.sequences = (char**)malloc(g_data.num_sequences * sizeof(char*));
    g_data.sequence_lengths = (int*)malloc(g_data.num_sequences * sizeof(int));

    const char alphabet[] = "ACGT";
    // Generate random DNA sequences
    for (int i = 0; i < g_data.num_sequences; ++i) {
        int len_variation = g_data.average_sequence_length / 4;
        int len = g_data.average_sequence_length + (mt_rand() % (2 * len_variation + 1)) - len_variation;
        if (len < g_data.motif_length) len = g_data.motif_length;
        g_data.sequence_lengths[i] = len;
        g_data.sequences[i] = (char*)malloc((len + 1) * sizeof(char));
        for (int j = 0; j < len; ++j) {
            g_data.sequences[i][j] = alphabet[mt_rand() % 4];
        }
        g_data.sequences[i][len] = '\0';
    }

    // Allocate other data structures
    g_data.motif_starts = (int*)malloc(g_data.num_sequences * sizeof(int));
    g_data.score_buffer = (double*)malloc(g_data.average_sequence_length * 2 * sizeof(double)); // Ample buffer size

    g_data.profile_matrix = (double**)malloc(g_data.motif_length * sizeof(double*));
    for (int i = 0; i < g_data.motif_length; ++i) {
        g_data.profile_matrix[i] = (double*)malloc(4 * sizeof(double));
    }

    g_data.final_result = 0;
}

void run_computation() {
    // 1. Initial random assignment of motif start sites
    for (int i = 0; i < g_data.num_sequences; ++i) {
        g_data.motif_starts[i] = mt_rand() % (g_data.sequence_lengths[i] - g_data.motif_length + 1);
    }

    // 2. Gibbs sampling main loop
    for (int iter = 0; iter < g_data.num_iterations; ++iter) {
        int seq_to_remove = iter % g_data.num_sequences;

        // a. Build profile matrix from all sequences EXCEPT seq_to_remove
        for (int i = 0; i < g_data.motif_length; ++i) {
            for (int j = 0; j < 4; ++j) {
                g_data.profile_matrix[i][j] = 1.0; // Pseudocount
            }
        }

        for (int i = 0; i < g_data.num_sequences; ++i) {
            if (i == seq_to_remove) continue;
            int start = g_data.motif_starts[i];
            for (int j = 0; j < g_data.motif_length; ++j) {
                char nucleotide = g_data.sequences[i][start + j];
                g_data.profile_matrix[j][nucleotide_to_index(nucleotide)]++;
            }
        }

        // Convert counts to probabilities
        for (int i = 0; i < g_data.motif_length; ++i) {
            double row_sum = 0.0;
            for (int j = 0; j < 4; ++j) {
                row_sum += g_data.profile_matrix[i][j];
            }
            for (int j = 0; j < 4; ++j) {
                g_data.profile_matrix[i][j] /= row_sum;
            }
        }

        // b. Score all possible motif positions in the removed sequence
        double total_score = 0.0;
        int max_start_pos = g_data.sequence_lengths[seq_to_remove] - g_data.motif_length;
        for (int s = 0; s <= max_start_pos; ++s) {
            double current_prob = 1.0;
            for (int j = 0; j < g_data.motif_length; ++j) {
                char nucleotide = g_data.sequences[seq_to_remove][s + j];
                current_prob *= g_data.profile_matrix[j][nucleotide_to_index(nucleotide)];
            }
            g_data.score_buffer[s] = current_prob;
            total_score += current_prob;
        }

        // c. Sample a new starting position based on weighted probabilities
        int new_start = 0;
        if (total_score < 1e-9) { // Failsafe for zero probability case
            new_start = mt_rand() % (max_start_pos + 1);
        } else {
            double rand_val = ((double)mt_rand() / (double)UINT32_MAX) * total_score;
            double cumulative_prob = 0.0;
            for (int s = 0; s <= max_start_pos; ++s) {
                cumulative_prob += g_data.score_buffer[s];
                if (rand_val <= cumulative_prob) {
                    new_start = s;
                    break;
                }
            }
             // Ensure new_start is valid due to potential float inaccuracies
            if (new_start > max_start_pos) new_start = max_start_pos;
        }

        // d. Update the motif start for the sequence
        g_data.motif_starts[seq_to_remove] = new_start;
    }

    // 3. Calculate a final result to prevent dead code elimination
    long long result_sum = 0;
    for (int i = 0; i < g_data.num_sequences; ++i) {
        result_sum += g_data.motif_starts[i];
    }
    g_data.final_result = result_sum;
}

void cleanup() {
    for (int i = 0; i < g_data.num_sequences; ++i) {
        free(g_data.sequences[i]);
    }
    free(g_data.sequences);
    free(g_data.sequence_lengths);
    free(g_data.motif_starts);
    free(g_data.score_buffer);

    for (int i = 0; i < g_data.motif_length; ++i) {
        free(g_data.profile_matrix[i]);
    }
    free(g_data.profile_matrix);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout
    printf("%lld\n", g_data.final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
