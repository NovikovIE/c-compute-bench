#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>

// --- MERSENNE TWISTER (DO NOT MODIFY) ---
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

// --- BENCHMARK DATA & PARAMETERS ---

typedef struct PatternNode {
    int* pattern;
    int length;
    struct PatternNode* next;
} PatternNode;

int num_sequences;
int avg_sequence_length;
int num_unique_items;
float min_support_threshold;

int** sequences_db;
int* sequence_lengths;
long long final_result;

// --- HELPER FUNCTIONS ---

void free_pattern_list(PatternNode* head) {
    PatternNode* current = head;
    while (current != NULL) {
        PatternNode* next = current->next;
        free(current->pattern);
        free(current);
        current = next;
    }
}

bool is_subsequence(const int* pattern, int pattern_len, const int* sequence, int sequence_len) {
    if (pattern_len == 0) return true;
    if (pattern_len > sequence_len) return false;

    int pattern_idx = 0;
    for (int seq_idx = 0; seq_idx < sequence_len; ++seq_idx) {
        if (sequence[seq_idx] == pattern[pattern_idx]) {
            pattern_idx++;
            if (pattern_idx == pattern_len) {
                return true;
            }
        }
    }
    return false;
}

// --- BENCHMARK IMPLEMENTATION ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s num_sequences avg_sequence_length num_unique_items min_support_threshold seed\n", argv[0]);
        exit(1);
    }

    num_sequences = atoi(argv[1]);
    avg_sequence_length = atoi(argv[2]);
    num_unique_items = atoi(argv[3]);
    min_support_threshold = atof(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    mt_seed(seed);

    sequences_db = (int**)malloc(num_sequences * sizeof(int*));
    sequence_lengths = (int*)malloc(num_sequences * sizeof(int));

    for (int i = 0; i < num_sequences; ++i) {
        int length = avg_sequence_length / 2 + (mt_rand() % avg_sequence_length);
        if (length <= 0) length = 1;
        sequence_lengths[i] = length;
        sequences_db[i] = (int*)malloc(length * sizeof(int));
        for (int j = 0; j < length; ++j) {
            sequences_db[i][j] = mt_rand() % num_unique_items;
        }
    }

    final_result = 0;
}

void run_computation() {
    int min_support_count = (int)(min_support_threshold * num_sequences);
    if (min_support_count < 1) min_support_count = 1;

    // Phase 1: Find frequent 1-item sequences (L1)
    PatternNode* frequent_L1 = NULL;
    int* item_counts = (int*)calloc(num_unique_items, sizeof(int));
    for (int i = 0; i < num_sequences; ++i) {
        // To avoid overcounting, use a flag for each item in a sequence
        bool* seen_in_sequence = (bool*)calloc(num_unique_items, sizeof(bool));
        for (int j = 0; j < sequence_lengths[i]; ++j) {
            int item = sequences_db[i][j];
            if (!seen_in_sequence[item]) {
                item_counts[item]++;
                seen_in_sequence[item] = true;
            }
        }
        free(seen_in_sequence);
    }

    for (int i = 0; i < num_unique_items; ++i) {
        if (item_counts[i] >= min_support_count) {
            PatternNode* node = (PatternNode*)malloc(sizeof(PatternNode));
            node->length = 1;
            node->pattern = (int*)malloc(sizeof(int));
            node->pattern[0] = i;
            node->next = frequent_L1;
            frequent_L1 = node;
        }
    }
    free(item_counts);

    if (frequent_L1 == NULL) {
        final_result = 0;
        return;
    }

    PatternNode* L_km1 = frequent_L1;
    int k = 2;

    while (1) {
        // Candidate Generation (Ck from L_{k-1} and L1)
        PatternNode* C_k = NULL;
        for (PatternNode* p_node = L_km1; p_node != NULL; p_node = p_node->next) {
            for (PatternNode* i_node = frequent_L1; i_node != NULL; i_node = i_node->next) {
                if (i_node->pattern[0] > p_node->pattern[p_node->length - 1]) {
                    PatternNode* candidate = (PatternNode*)malloc(sizeof(PatternNode));
                    candidate->length = k;
                    candidate->pattern = (int*)malloc(k * sizeof(int));
                    memcpy(candidate->pattern, p_node->pattern, p_node->length * sizeof(int));
                    candidate->pattern[k - 1] = i_node->pattern[0];
                    candidate->next = C_k;
                    C_k = candidate;
                }
            }
        }
        
        // Support Counting & Pruning (Generate Lk from Ck)
        PatternNode* L_k = NULL;
        PatternNode* current_candidate = C_k;
        while (current_candidate != NULL) {
            int support = 0;
            for (int i = 0; i < num_sequences; ++i) {
                if (is_subsequence(current_candidate->pattern, current_candidate->length, sequences_db[i], sequence_lengths[i])) {
                    support++;
                }
            }

            if (support >= min_support_count) {
                PatternNode* frequent_node = (PatternNode*)malloc(sizeof(PatternNode));
                frequent_node->length = current_candidate->length;
                frequent_node->pattern = (int*)malloc(frequent_node->length * sizeof(int));
                memcpy(frequent_node->pattern, current_candidate->pattern, frequent_node->length * sizeof(int));
                frequent_node->next = L_k;
                L_k = frequent_node;
            }
            current_candidate = current_candidate->next;
        }
        free_pattern_list(C_k);

        if (L_k == NULL) {
            // No more frequent patterns, L_km1 is the final set
            break;
        }

        free_pattern_list(L_km1);
        L_km1 = L_k;
        k++;
    }

    // Accumulate result from the final pattern set (L_km1)
    final_result = 0;
    for (PatternNode* node = L_km1; node != NULL; node = node->next) {
        for (int i = 0; i < node->length; ++i) {
            final_result += node->pattern[i];
        }
    }

    free_pattern_list(L_km1);
}

void cleanup() {
    for (int i = 0; i < num_sequences; ++i) {
        free(sequences_db[i]);
    }
    free(sequences_db);
    free(sequence_lengths);
}

// --- MAIN FUNCTION ---
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