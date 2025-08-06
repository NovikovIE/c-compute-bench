#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

// --- Mersenne Twister (MT19937) ---
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

// --- Data Structures ---

typedef struct {
    int* items;
    int length;
} Sequence;

typedef struct {
    Sequence* sequences;
    int num_sequences;
    int avg_sequence_length;
    int num_unique_items;
    int min_support;
} Database;

typedef struct {
    int seq_idx; 
    int pos;     
} ProjectedSequence;

// --- Global Variables ---
Database* g_db;
long long g_frequent_pattern_count;

// --- Forward Declarations ---
void recursive_mine(ProjectedSequence* projected_db, int db_size);

// --- Benchmark Functions ---

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <num_sequences> <avg_sequence_length> <num_unique_items> <min_support_threshold> <seed>\n", argv[0]);
        exit(1);
    }

    g_db = (Database*)malloc(sizeof(Database));
    if (!g_db) {
        perror("Failed to allocate memory for Database");
        exit(1);
    }

    g_db->num_sequences = atoi(argv[1]);
    g_db->avg_sequence_length = atoi(argv[2]);
    g_db->num_unique_items = atoi(argv[3]);
    g_db->min_support = atoi(argv[4]);
    uint32_t seed = atoi(argv[5]);

    mt_seed(seed);

    g_db->sequences = (Sequence*)malloc(g_db->num_sequences * sizeof(Sequence));
    if (!g_db->sequences) {
        perror("Failed to allocate memory for sequences array");
        exit(1);
    }

    for (int i = 0; i < g_db->num_sequences; ++i) {
        int length = g_db->avg_sequence_length > 1 ? (g_db->avg_sequence_length / 2) + (mt_rand() % g_db->avg_sequence_length) : 1;
        g_db->sequences[i].length = length;
        g_db->sequences[i].items = (int*)malloc(length * sizeof(int));
        if (!g_db->sequences[i].items) {
            perror("Failed to allocate memory for sequence items");
            exit(1);
        }
        for (int j = 0; j < length; ++j) {
            g_db->sequences[i].items[j] = mt_rand() % g_db->num_unique_items;
        }
    }
}

void run_computation() {
    g_frequent_pattern_count = 0;

    int* item_counts = (int*)calloc(g_db->num_unique_items, sizeof(int));
    int* last_seen_in_seq = (int*)malloc(g_db->num_unique_items * sizeof(int));
    if (!item_counts || !last_seen_in_seq) { exit(1); }
    
    // Step 1: Find frequent items of length 1 (1-patterns)
    for (int i = 0; i < g_db->num_sequences; i++) {
        // Reset seen items for this new sequence
        memset(last_seen_in_seq, -1, g_db->num_unique_items * sizeof(int));
        for (int j = 0; j < g_db->sequences[i].length; j++) {
            int item = g_db->sequences[i].items[j];
            if (last_seen_in_seq[item] != i) {
                item_counts[item]++;
                last_seen_in_seq[item] = i; // Mark item as seen in this sequence
            }
        }
    }
    
    // Step 2: For each frequent 1-pattern, create a projected DB and recurse
    for (int item = 0; item < g_db->num_unique_items; item++) {
        if (item_counts[item] >= g_db->min_support) {
            g_frequent_pattern_count++;

            ProjectedSequence* projected_db = (ProjectedSequence*)malloc(item_counts[item] * sizeof(ProjectedSequence));
            if (!projected_db) continue;
            
            int p_db_idx = 0;
            for (int i = 0; i < g_db->num_sequences; i++) {
                for (int j = 0; j < g_db->sequences[i].length; j++) {
                    if (g_db->sequences[i].items[j] == item) {
                        projected_db[p_db_idx++] = (ProjectedSequence){.seq_idx = i, .pos = j + 1};
                        break; // Move to the next sequence
                    }
                }
            }
            recursive_mine(projected_db, p_db_idx);
            free(projected_db);
        }
    }

    free(item_counts);
    free(last_seen_in_seq);
}

void recursive_mine(ProjectedSequence* projected_db, int db_size) {
    if (db_size < g_db->min_support) return;

    int* local_item_counts = (int*)calloc(g_db->num_unique_items, sizeof(int));
    if (!local_item_counts) exit(1);

    // To count suffix items once per projected sequence
    int* temp_last_seen = (int*)calloc(g_db->num_unique_items, sizeof(int));
    if (!temp_last_seen) { free(local_item_counts); exit(1); }

    // Pass 1: Count frequencies of items following the current pattern
    for (int i = 0; i < db_size; i++) {
        int seq_idx = projected_db[i].seq_idx;
        int start_pos = projected_db[i].pos;
        for (int j = start_pos; j < g_db->sequences[seq_idx].length; j++) {
            int item = g_db->sequences[seq_idx].items[j];
            if (temp_last_seen[item] != i + 1) { // Use i+1 to distinguish from 0
                local_item_counts[item]++;
                temp_last_seen[item] = i + 1;
            }
        }
    }

    // Pass 2: For each locally frequent item, build its projected DB and recurse
    for (int item = 0; item < g_db->num_unique_items; item++) {
        if (local_item_counts[item] >= g_db->min_support) {
            g_frequent_pattern_count++;

            ProjectedSequence* next_projected_db = (ProjectedSequence*)malloc(local_item_counts[item] * sizeof(ProjectedSequence));
            if (!next_projected_db)  continue;

            int next_db_idx = 0;
            for (int i = 0; i < db_size; i++) {
                int seq_idx = projected_db[i].seq_idx;
                int start_pos = projected_db[i].pos;
                for (int j = start_pos; j < g_db->sequences[seq_idx].length; j++) {
                    if (g_db->sequences[seq_idx].items[j] == item) {
                        next_projected_db[next_db_idx++] = (ProjectedSequence){.seq_idx = seq_idx, .pos = j + 1};
                        break;
                    }
                }
            }
            recursive_mine(next_projected_db, next_db_idx);
            free(next_projected_db);
        }
    }
    
    free(local_item_counts);
    free(temp_last_seen);
}


void cleanup() {
    for (int i = 0; i < g_db->num_sequences; ++i) {
        free(g_db->sequences[i].items);
    }
    free(g_db->sequences);
    free(g_db);
}

int main(int argc, char* argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%lld\n", g_frequent_pattern_count);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
