#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- BEGIN MERSENNE TWISTER ---
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

// --- BENCHMARK DATA AND PARAMETERS ---
typedef struct {
    int* items;
    int size;
    float support;
} Itemset;

// Parameters
static int P_NUM_FREQUENT_ITEMSETS;
static float P_MIN_CONFIDENCE_THRESHOLD;

// Data structures
static Itemset* frequent_itemsets;
static unsigned int final_rule_count = 0;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_frequent_itemsets> <min_confidence_threshold> <seed>\n", argv[0]);
        exit(1);
    }

    P_NUM_FREQUENT_ITEMSETS = atoi(argv[1]);
    P_MIN_CONFIDENCE_THRESHOLD = atof(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);
    mt_seed(seed);
    
    frequent_itemsets = (Itemset*)malloc(P_NUM_FREQUENT_ITEMSETS * sizeof(Itemset));
    if (!frequent_itemsets) {
        fprintf(stderr, "Failed to allocate memory for frequent_itemsets\n");
        exit(1);
    }

    for (int i = 0; i < P_NUM_FREQUENT_ITEMSETS; i++) {
        // Generate itemsets of size k, where 10 <= k <= 18.
        // This range provides a computationally intensive but tractable number of subsets (2^k - 2).
        int k = 10 + (mt_rand() % 9);
        frequent_itemsets[i].size = k;
        
        // Frequent itemsets have relatively low support in practice.
        frequent_itemsets[i].support = ((float)mt_rand() / (float)UINT32_MAX) * 0.2f + 0.05f; // Support in [0.05, 0.25]
        
        frequent_itemsets[i].items = (int*)malloc(k * sizeof(int));
        if (!frequent_itemsets[i].items) {
            fprintf(stderr, "Failed to allocate memory for an itemset\n");
            exit(1);
        }

        // Generate a sorted list of unique items to represent the itemset.
        // Start with a random base and add random increments.
        int current_item_id = mt_rand() % 1000;
        for (int j = 0; j < k; j++) {
            frequent_itemsets[i].items[j] = current_item_id;
            current_item_id += (mt_rand() % 10) + 1;
        }
    }
}

void run_computation() {
    unsigned int local_rule_count = 0;
    for (int i = 0; i < P_NUM_FREQUENT_ITEMSETS; i++) {
        Itemset current_set = frequent_itemsets[i];
        int k = current_set.size;

        // Total number of non-empty proper subsets is 2^k - 2.
        // We use a bitmask to iterate through all of them.
        unsigned long num_subsets = 1UL << k;
        
        // Loop from 1 to 2^k - 2 to represent all antecedents.
        for (unsigned long j = 1; j < num_subsets - 1; j++) {
            // In a real system, we'd look up the support of the antecedent subset.
            // Here, we simulate this lookup with a random number generation.
            // The support of a subset (antecedent) must be >= support of the superset.
            float antecedent_support = current_set.support + 
                                     ((float)mt_rand() / (float)UINT32_MAX) * (1.0f - current_set.support);

            // Calculate confidence: Confidence(A -> B) = Support(A U B) / Support(A)
            float confidence = current_set.support / antecedent_support;

            if (confidence >= P_MIN_CONFIDENCE_THRESHOLD) {
                local_rule_count++;
            }
        }
    }
    final_rule_count = local_rule_count;
}

void cleanup() {
    for (int i = 0; i < P_NUM_FREQUENT_ITEMSETS; i++) {
        free(frequent_itemsets[i].items);
    }
    free(frequent_itemsets);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    cleanup();

    printf("%u\n", final_rule_count);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
