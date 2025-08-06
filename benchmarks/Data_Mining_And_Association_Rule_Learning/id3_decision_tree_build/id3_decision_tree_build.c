#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <float.h>

// --- Mersenne Twister (MT19937) Generator - DO NOT MODIFY ---
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

// --- Benchmark Globals ---
int NUM_SAMPLES;
int NUM_FEATURES;
int NUM_CLASSES;
int NUM_FEATURE_VALUES;

// Dataset
int** data; // Shape: [NUM_SAMPLES][NUM_FEATURES]
int* labels; // Shape: [NUM_SAMPLES]

// To track used features during tree construction simulation
int* feature_used;

// Final result to prevent dead-code elimination
int final_result;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s num_samples num_features num_classes num_feature_values seed\n", argv[0]);
        exit(1);
    }

    NUM_SAMPLES = atoi(argv[1]);
    NUM_FEATURES = atoi(argv[2]);
    NUM_CLASSES = atoi(argv[3]);
    NUM_FEATURE_VALUES = atoi(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    mt_seed(seed);

    // Allocate data matrix
    data = (int**)malloc(NUM_SAMPLES * sizeof(int*));
    if (!data) { perror("malloc failed"); exit(1); }
    for (int i = 0; i < NUM_SAMPLES; ++i) {
        data[i] = (int*)malloc(NUM_FEATURES * sizeof(int));
        if (!data[i]) { perror("malloc failed"); exit(1); }
    }

    // Allocate labels array
    labels = (int*)malloc(NUM_SAMPLES * sizeof(int));
    if (!labels) { perror("malloc failed"); exit(1); }

    // Allocate feature usage tracker
    feature_used = (int*)calloc(NUM_FEATURES, sizeof(int));
    if (!feature_used) { perror("malloc failed"); exit(1); }

    // Populate with random data
    for (int i = 0; i < NUM_SAMPLES; ++i) {
        for (int j = 0; j < NUM_FEATURES; ++j) {
            data[i][j] = mt_rand() % NUM_FEATURE_VALUES;
        }
        labels[i] = mt_rand() % NUM_CLASSES;
    }

    final_result = 0;
}

void run_computation() {
    // Simulate building a decision tree by iteratively selecting the best feature
    // to split on, up to a fixed depth. This captures the core computation of ID3.
    const int MAX_TREE_DEPTH = 10;

    // Pre-calculate the total entropy of the entire dataset (H(S))
    double total_entropy = 0.0;
    int* total_class_counts = (int*)calloc(NUM_CLASSES, sizeof(int));
    if (!total_class_counts) { perror("calloc failed"); exit(1); }
    for (int i = 0; i < NUM_SAMPLES; i++) {
        total_class_counts[labels[i]]++;
    }
    for (int i = 0; i < NUM_CLASSES; i++) {
        if (total_class_counts[i] > 0) {
            double p = (double)total_class_counts[i] / NUM_SAMPLES;
            total_entropy -= p * log2(p);
        }
    }
    free(total_class_counts);

    // Main loop to simulate finding the best split for each level of the tree
    for (int depth = 0; depth < MAX_TREE_DEPTH; depth++) {
        double max_gain = -1.0;
        int best_feature = -1;

        // Find the best feature among the unused ones
        for (int f = 0; f < NUM_FEATURES; f++) {
            if (feature_used[f]) {
                continue;
            }
            
            // Calculate Information Gain for feature 'f'.
            // This requires partitioning data by 'f' and calculating the weighted entropy.
            int** partition_class_counts = (int**)malloc(NUM_FEATURE_VALUES * sizeof(int*));
            for(int i = 0; i < NUM_FEATURE_VALUES; ++i) {
                partition_class_counts[i] = (int*)calloc(NUM_CLASSES, sizeof(int));
            }
            int* partition_total_counts = (int*)calloc(NUM_FEATURE_VALUES, sizeof(int));

            // Create distributions for each partition
            for (int s = 0; s < NUM_SAMPLES; s++) {
                int value = data[s][f];
                int class_label = labels[s];
                partition_class_counts[value][class_label]++;
                partition_total_counts[value]++;
            }

            // Calculate the weighted average entropy of the partitions
            double weighted_entropy = 0.0;
            for (int v = 0; v < NUM_FEATURE_VALUES; v++) {
                if (partition_total_counts[v] == 0) continue;

                double partition_entropy = 0.0;
                for (int c = 0; c < NUM_CLASSES; c++) {
                    if (partition_class_counts[v][c] > 0) {
                        double p = (double)partition_class_counts[v][c] / partition_total_counts[v];
                        partition_entropy -= p * log2(p);
                    }
                }
                weighted_entropy += ((double)partition_total_counts[v] / NUM_SAMPLES) * partition_entropy;
            }

            double gain = total_entropy - weighted_entropy;

            if (gain > max_gain) {
                max_gain = gain;
                best_feature = f;
            }
            
            // Cleanup for this feature's calculation
            for(int i = 0; i < NUM_FEATURE_VALUES; ++i) {
                free(partition_class_counts[i]);
            }
            free(partition_class_counts);
            free(partition_total_counts);
        }

        if (best_feature != -1) {
            feature_used[best_feature] = 1;
            final_result += best_feature; // Accumulate result
        } else {
            break; // No more features with positive gain
        }
    }
}

void cleanup() {
    for (int i = 0; i < NUM_SAMPLES; ++i) {
        free(data[i]);
    }
    free(data);
    free(labels);
    free(feature_used);
}

int main(int argc, char* argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout
    printf("%d\n", final_result);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
