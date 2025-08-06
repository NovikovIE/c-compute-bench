#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

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

// Benchmark parameters
static int NUM_TRANSACTIONS;
static int NUM_UNIQUE_ITEMS;
static int MIN_SUPPORT_THRESHOLD;

// Data structures: Transactions are represented as bitmasks
static uint64_t *transactions;

// Result
static int total_frequent_itemsets_result;

// Helper for qsort
int compare_u64(const void *a, const void *b) {
    uint64_t val_a = *(const uint64_t *)a;
    uint64_t val_b = *(const uint64_t *)b;
    if (val_a < val_b) return -1;
    if (val_a > val_b) return 1;
    return 0;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_transactions num_unique_items min_support_threshold seed\n", argv[0]);
        exit(1);
    }

    NUM_TRANSACTIONS = atoi(argv[1]);
    NUM_UNIQUE_ITEMS = atoi(argv[2]);
    MIN_SUPPORT_THRESHOLD = atoi(argv[3]);
    uint32_t seed = atoi(argv[4]);
    
    if (NUM_UNIQUE_ITEMS > 64) {
        fprintf(stderr, "Error: num_unique_items cannot exceed 64 for this bitmask implementation.\n");
        exit(1);
    }

    mt_seed(seed);

    transactions = (uint64_t *)malloc(NUM_TRANSACTIONS * sizeof(uint64_t));
    if (!transactions) {
        perror("Failed to allocate memory for transactions");
        exit(1);
    }

    for (int i = 0; i < NUM_TRANSACTIONS; ++i) {
        transactions[i] = 0;
        int items_in_basket = 3 + (mt_rand() % 5); // 3 to 7 items
        for (int j = 0; j < items_in_basket; ++j) {
            // Generate items with a skewed distribution (low-ID items are more common)
            double r = (double)mt_rand() / (double)UINT32_MAX;
            int item_id = (int)((double)NUM_UNIQUE_ITEMS * (1.0 - sqrt(1.0 - r)));
            if (item_id >= NUM_UNIQUE_ITEMS) item_id = NUM_UNIQUE_ITEMS - 1;
            transactions[i] |= (1ULL << item_id);
        }
    }
}

void run_computation() {
    total_frequent_itemsets_result = 0;

    // --- Level 1: Find frequent 1-itemsets (L1) ---
    int *item_counts = (int *)calloc(NUM_UNIQUE_ITEMS, sizeof(int));
    if (!item_counts) { perror("calloc"); exit(1); }

    for (int i = 0; i < NUM_TRANSACTIONS; ++i) {
        for (int j = 0; j < NUM_UNIQUE_ITEMS; ++j) {
            if ((transactions[i] >> j) & 1) {
                item_counts[j]++;
            }
        }
    }

    uint64_t *L_current = (uint64_t *)malloc(NUM_UNIQUE_ITEMS * sizeof(uint64_t));
    if (!L_current) { perror("malloc"); exit(1); }
    int L_current_size = 0;
    for (int i = 0; i < NUM_UNIQUE_ITEMS; ++i) {
        if (item_counts[i] >= MIN_SUPPORT_THRESHOLD) {
            L_current[L_current_size++] = (1ULL << i);
        }
    }
    free(item_counts);
    total_frequent_itemsets_result += L_current_size;

    // --- Main Loop: Generate L(k) from L(k-1) ---
    for (int k = 2; L_current_size > 0 && k <= NUM_UNIQUE_ITEMS; ++k) {
        // --- Generate Candidates Ck ---
        size_t max_candidates = (L_current_size * (L_current_size - 1)) / 2;
        if (max_candidates == 0) break;
        uint64_t *candidates = (uint64_t *)malloc(max_candidates * sizeof(uint64_t));
        if (!candidates) { perror("malloc"); exit(1); }
        int Ck_size = 0;

        for (int i = 0; i < L_current_size; ++i) {
            for (int j = i + 1; j < L_current_size; ++j) {
                uint64_t candidate = L_current[i] | L_current[j];
                // Using popcount to ensure we are joining k-1 itemsets to form a k-itemset
                // This is a simplified join step, not the classic Apriori prefix-based join
                if (__builtin_popcountll(candidate) == k) {
                    candidates[Ck_size++] = candidate;
                }
            }
        }

        // Sort and remove duplicates to get the final candidate set Ck
        qsort(candidates, Ck_size, sizeof(uint64_t), compare_u64);
        int unique_Ck_size = 0;
        if (Ck_size > 0) {
            unique_Ck_size = 1;
            for (int i = 1; i < Ck_size; ++i) {
                if (candidates[i] != candidates[i - 1]) {
                    candidates[unique_Ck_size++] = candidates[i];
                }
            }
        }
        Ck_size = unique_Ck_size;

        // --- Count Support for Candidates Ck ---
        int *Ck_counts = (int *)calloc(Ck_size, sizeof(int));
        if (!Ck_counts) { perror("calloc"); exit(1); }

        for (int i = 0; i < NUM_TRANSACTIONS; ++i) {
            for (int j = 0; j < Ck_size; ++j) {
                if ((transactions[i] & candidates[j]) == candidates[j]) {
                    Ck_counts[j]++;
                }
            }
        }

        // --- Filter Ck to get Lk ---
        uint64_t *L_next = (uint64_t *)malloc(Ck_size * sizeof(uint64_t));
        if (!L_next) { perror("malloc"); exit(1); }
        int L_next_size = 0;
        for (int i = 0; i < Ck_size; ++i) {
            if (Ck_counts[i] >= MIN_SUPPORT_THRESHOLD) {
                L_next[L_next_size++] = candidates[i];
            }
        }

        // Cleanup for this level and prepare for the next
        free(candidates);
        free(Ck_counts);
        free(L_current);

        L_current = L_next;
        L_current_size = L_next_size;
        total_frequent_itemsets_result += L_current_size;
    }

    free(L_current); // Free the final level's list
}

void cleanup() {
    if (transactions) {
        free(transactions);
    }
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
    printf("%d\n", total_frequent_itemsets_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
