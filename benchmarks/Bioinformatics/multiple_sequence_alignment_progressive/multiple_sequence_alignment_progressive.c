#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>

// --- Mersenne Twister (verbatim) ---
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
// --- End Mersenne Twister ---

// --- Benchmark Globals and Defs ---
#define MAX(a,b) (((a)>(b))?(a):(b))

const int MATCH_SCORE = 5;
const int MISMATCH_PENALTY = -4;
const int GAP_PENALTY = -6;

typedef struct {
    int num_sequences;
    int avg_sequence_length;
    char** sequences;
    int* sequence_lengths;
    int max_len;
    long long final_score;
} BenchmarkData;

static BenchmarkData g_data;

typedef struct {
    int index;
    int score;
} IndexScore;

// --- Function Prototypes ---
void setup_benchmark(int argc, char *argv[]);
void run_computation();
void cleanup();
int needleman_wunsch_score(const char* seq1, int len1, const char* seq2, int len2, int** dp_table);
int compare_index_score(const void* a, const void* b);

// --- Benchmark Implementation ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_sequences> <average_sequence_length> <seed>\n", argv[0]);
        exit(1);
    }
    g_data.num_sequences = atoi(argv[1]);
    g_data.avg_sequence_length = atoi(argv[2]);
    uint32_t seed = (uint32_t)strtoul(argv[3], NULL, 10);
    mt_seed(seed);

    if (g_data.num_sequences <= 1 || g_data.avg_sequence_length <= 0) {
        fprintf(stderr, "FATAL: Invalid parameters.\n");
        exit(1);
    }

    g_data.sequences = (char**)malloc(g_data.num_sequences * sizeof(char*));
    g_data.sequence_lengths = (int*)malloc(g_data.num_sequences * sizeof(int));
    if (!g_data.sequences || !g_data.sequence_lengths) {
        fprintf(stderr, "FATAL: Memory allocation failed for sequence pointers.\n");
        exit(1);
    }
    
    g_data.max_len = 0;
    const char alphabet[] = "ACGT";
    const int alphabet_size = 4;

    for (int i = 0; i < g_data.num_sequences; i++) {
        int len = (g_data.avg_sequence_length / 2) + (mt_rand() % g_data.avg_sequence_length);
        if (len == 0) len = 1;
        g_data.sequence_lengths[i] = len;
        if (len > g_data.max_len) {
            g_data.max_len = len;
        }

        g_data.sequences[i] = (char*)malloc((len + 1) * sizeof(char));
        if (!g_data.sequences[i]) {
            fprintf(stderr, "FATAL: Memory allocation failed for a sequence.\n");
            exit(1);
        }

        for (int j = 0; j < len; j++) {
            g_data.sequences[i][j] = alphabet[mt_rand() % alphabet_size];
        }
        g_data.sequences[i][len] = '\0';
    }
}

int needleman_wunsch_score(const char* seq1, int len1, const char* seq2, int len2, int** dp_table) {
    for (int i = 0; i <= len1; i++) {
        dp_table[i][0] = i * GAP_PENALTY;
    }
    for (int j = 0; j <= len2; j++) {
        dp_table[0][j] = j * GAP_PENALTY;
    }

    for (int i = 1; i <= len1; i++) {
        for (int j = 1; j <= len2; j++) {
            int score_diag = dp_table[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? MATCH_SCORE : MISMATCH_PENALTY);
            int score_up = dp_table[i - 1][j] + GAP_PENALTY;
            int score_left = dp_table[i][j - 1] + GAP_PENALTY;
            dp_table[i][j] = MAX(score_diag, MAX(score_up, score_left));
        }
    }
    return dp_table[len1][len2];
}

int compare_index_score(const void* a, const void* b) {
    const IndexScore* is_a = (const IndexScore*)a;
    const IndexScore* is_b = (const IndexScore*)b;
    if (is_b->score > is_a->score) return 1;
    if (is_b->score < is_a->score) return -1;
    return 0;
}

void run_computation() {
    int n = g_data.num_sequences;
    
    int dp_dim = g_data.max_len + 1;
    int** dp_table = (int**)malloc(dp_dim * sizeof(int*));
    for (int i = 0; i < dp_dim; i++) {
        dp_table[i] = (int*)malloc(dp_dim * sizeof(int));
    }

    int** scores = (int**)malloc(n * sizeof(int*));
    for (int i = 0; i < n; i++) {
        scores[i] = (int*)malloc(n * sizeof(int));
    }

    for (int i = 0; i < n; i++) {
        scores[i][i] = 0;
        for (int j = i + 1; j < n; j++) {
            int score = needleman_wunsch_score(g_data.sequences[i], g_data.sequence_lengths[i],
                                               g_data.sequences[j], g_data.sequence_lengths[j],
                                               dp_table);
            scores[i][j] = scores[j][i] = score;
        }
    }

    int center_idx = 0;
    long max_total_score = -1000000000000LL; 
    for (int i = 0; i < n; i++) {
        long current_total_score = 0;
        for (int j = 0; j < n; j++) {
            current_total_score += scores[i][j];
        }
        if (current_total_score > max_total_score) {
            max_total_score = current_total_score;
            center_idx = i;
        }
    }

    IndexScore* order_to_align = (IndexScore*)malloc((n - 1) * sizeof(IndexScore));
    int k = 0;
    for (int i = 0; i < n; i++) {
        if (i == center_idx) continue;
        order_to_align[k].index = i;
        order_to_align[k].score = scores[center_idx][i];
        k++;
    }
    qsort(order_to_align, n - 1, sizeof(IndexScore), compare_index_score);

    long long total_alignment_score = 0;
    int* profile_members = (int*)malloc(n * sizeof(int));
    int profile_size = 0;
    
    profile_members[profile_size++] = center_idx;

    for (int i = 0; i < n - 1; i++) {
        int new_seq_idx = order_to_align[i].index;
        long long current_alignment_score = 0;
        
        for (int p = 0; p < profile_size; p++) {
            int profile_seq_idx = profile_members[p];
            current_alignment_score += scores[new_seq_idx][profile_seq_idx];
        }
        
        total_alignment_score += current_alignment_score / profile_size;
        
        profile_members[profile_size++] = new_seq_idx;
    }
    
    g_data.final_score = total_alignment_score;

    for (int i = 0; i < dp_dim; i++) {
        free(dp_table[i]);
    }
    free(dp_table);
    for (int i = 0; i < n; i++) {
        free(scores[i]);
    }
    free(scores);
    free(order_to_align);
    free(profile_members);
}

void cleanup() {
    for (int i = 0; i < g_data.num_sequences; i++) {
        free(g_data.sequences[i]);
    }
    free(g_data.sequences);
    free(g_data.sequence_lengths);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);
    
    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%lld\n", g_data.final_score);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
