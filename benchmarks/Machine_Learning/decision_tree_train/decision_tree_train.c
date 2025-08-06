#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <stdbool.h>
#include <float.h>
#include <string.h>

// --- MERSENNE TWISTER (MT19937) --- (DO NOT MODIFY)
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

// --- BENCHMARK SPECIFIC CODE ---
#define NUM_THRESHOLDS_TO_CHECK 10
#define MIN_SAMPLES_LEAF 5

// Data structures for the decision tree
typedef struct {
    int feature_index;    // Index of the feature to split on
    float threshold;      // Threshold value for the split
    int left_child_idx;   // Index of the left child node in the tree array
    int right_child_idx;  // Index of the right child node in the tree array
    int class_label;      // Predicted class label (if a leaf node)
    bool is_leaf;
} Node;

typedef struct {
    int num_samples;
    int num_features;
    int max_depth;
    float **features;     // num_samples x num_features
    int *labels;          // num_samples
    int *sample_indices;  // Initial indices for root node
} Dataset;

// Global state for the benchmark
struct {
    Dataset data;
    Node *tree;
    int node_count;
    int max_nodes;
} benchmark_state;

long final_result = 0;

// Function prototypes
int build_tree(int* current_indices, int num_node_samples, int depth);

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_samples num_features max_depth seed\n", argv[0]);
        exit(1);
    }

    benchmark_state.data.num_samples = atoi(argv[1]);
    benchmark_state.data.num_features = atoi(argv[2]);
    benchmark_state.data.max_depth = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    mt_seed(seed);

    // Allocate dataset
    benchmark_state.data.features = (float **)malloc(benchmark_state.data.num_samples * sizeof(float *));
    if (!benchmark_state.data.features) exit(1);
    for (int i = 0; i < benchmark_state.data.num_samples; i++) {
        benchmark_state.data.features[i] = (float *)malloc(benchmark_state.data.num_features * sizeof(float));
        if (!benchmark_state.data.features[i]) exit(1);
    }
    benchmark_state.data.labels = (int *)malloc(benchmark_state.data.num_samples * sizeof(int));
    if (!benchmark_state.data.labels) exit(1);
    benchmark_state.data.sample_indices = (int *)malloc(benchmark_state.data.num_samples * sizeof(int));
    if (!benchmark_state.data.sample_indices) exit(1);

    // Generate random data
    for (int i = 0; i < benchmark_state.data.num_samples; i++) {
        for (int j = 0; j < benchmark_state.data.num_features; j++) {
            benchmark_state.data.features[i][j] = (float)mt_rand() / (float)UINT32_MAX;
        }
        // Make label correlated with the first feature for a more realistic scenario
        benchmark_state.data.labels[i] = (benchmark_state.data.features[i][0] > 0.5f) ? 1 : 0;
        benchmark_state.data.sample_indices[i] = i;
    }
    
    // Allocate tree structure
    benchmark_state.max_nodes = (1 << (benchmark_state.data.max_depth + 1)) - 1;
    benchmark_state.tree = (Node *)malloc(benchmark_state.max_nodes * sizeof(Node));
    if (!benchmark_state.tree) exit(1);
    benchmark_state.node_count = 0;
}

// Calculate Gini impurity from class counts
float calculate_gini_from_counts(int count0, int count1) {
    int total = count0 + count1;
    if (total == 0) return 0.0f;
    float p0 = (float)count0 / total;
    float p1 = (float)count1 / total;
    return 1.0f - (p0 * p0 + p1 * p1);
}

// Recursively build the decision tree
int build_tree(int* current_indices, int num_node_samples, int depth) {
    if (benchmark_state.node_count >= benchmark_state.max_nodes) {
        // Fallback for safety, should not be hit with proper max_depth
        return -1;
    }

    int node_idx = benchmark_state.node_count++;
    Node *node = &benchmark_state.tree[node_idx];
    node->is_leaf = false;
    node->left_child_idx = -1;
    node->right_child_idx = -1;

    // Check for stopping conditions: pure node or max depth reached
    int class_counts[2] = {0, 0};
    for (int i = 0; i < num_node_samples; i++) {
        class_counts[benchmark_state.data.labels[current_indices[i]]]++;
    }

    if (depth == benchmark_state.data.max_depth || num_node_samples <= MIN_SAMPLES_LEAF || class_counts[0] == 0 || class_counts[1] == 0) {
        node->is_leaf = true;
        node->class_label = (class_counts[1] > class_counts[0]) ? 1 : 0;
        return node_idx;
    }

    // Find the best split (feature and threshold)
    float best_gini = 1.0f; 
    int best_feature = -1;
    float best_threshold = -1.0f;

    for (int f = 0; f < benchmark_state.data.num_features; f++) {
        for (int t = 0; t < NUM_THRESHOLDS_TO_CHECK; t++) {
            float threshold = (float)mt_rand() / (float)UINT32_MAX;
            int left_counts[2] = {0, 0}, right_counts[2] = {0, 0};

            for (int i = 0; i < num_node_samples; i++) {
                int sample_idx = current_indices[i];
                if (benchmark_state.data.features[sample_idx][f] < threshold) {
                    left_counts[benchmark_state.data.labels[sample_idx]]++;
                } else {
                    right_counts[benchmark_state.data.labels[sample_idx]]++;
                }
            }

            int num_left = left_counts[0] + left_counts[1];
            int num_right = right_counts[0] + right_counts[1];

            if (num_left == 0 || num_right == 0) continue;

            float gini_left = calculate_gini_from_counts(left_counts[0], left_counts[1]);
            float gini_right = calculate_gini_from_counts(right_counts[0], right_counts[1]);
            float weighted_gini = ((float)num_left / num_node_samples) * gini_left + 
                                  ((float)num_right / num_node_samples) * gini_right;

            if (weighted_gini < best_gini) {
                best_gini = weighted_gini;
                best_feature = f;
                best_threshold = threshold;
            }
        }
    }

    if (best_feature == -1) { // No split improved the impurity
        node->is_leaf = true;
        node->class_label = (class_counts[1] > class_counts[0]) ? 1 : 0;
        return node_idx;
    }

    node->feature_index = best_feature;
    node->threshold = best_threshold;

    // Partition samples for child nodes
    int *left_indices = (int*)malloc(num_node_samples * sizeof(int));
    int *right_indices = (int*)malloc(num_node_samples * sizeof(int));
    if (!left_indices || !right_indices) exit(1);

    int num_left = 0, num_right = 0;
    for (int i = 0; i < num_node_samples; i++) {
        int sample_idx = current_indices[i];
        if (benchmark_state.data.features[sample_idx][best_feature] < best_threshold) {
            left_indices[num_left++] = sample_idx;
        } else {
            right_indices[num_right++] = sample_idx;
        }
    }
    
    node->left_child_idx = build_tree(left_indices, num_left, depth + 1);
    node->right_child_idx = build_tree(right_indices, num_right, depth + 1);

    free(left_indices);
    free(right_indices);

    return node_idx;
}

void run_computation() {
    build_tree(benchmark_state.data.sample_indices, benchmark_state.data.num_samples, 0);

    // Accumulate a value from the tree to prevent dead-code elimination
    long hash = 0;
    for (int i = 0; i < benchmark_state.node_count; i++) {
        Node *n = &benchmark_state.tree[i];
        if (n->is_leaf) {
            hash = (hash * 31 + n->class_label);
        } else {
            hash = (hash * 37 + n->feature_index);
            hash = (hash * 41 + (long)(n->threshold * 100.0f)); 
        }
    }
    final_result = hash;
}

void cleanup() {
    for (int i = 0; i < benchmark_state.data.num_samples; i++) {
        free(benchmark_state.data.features[i]);
    }
    free(benchmark_state.data.features);
    free(benchmark_state.data.labels);
    free(benchmark_state.data.sample_indices);
    free(benchmark_state.tree);
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
    printf("%ld\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
