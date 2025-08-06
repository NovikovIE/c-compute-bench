#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <string.h>

// --- Mersenne Twister (MT19937) Generator ---
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

// --- Benchmark-specific Code ---

// --- Parameters ---
int P_NUM_SAMPLES;
int P_NUM_FEATURES;
int P_NUM_TREES;
int P_MAX_DEPTH;
int P_FEATURES_PER_SPLIT;

// --- Data Structures ---
typedef struct {
    int feature_index;
    float threshold;
    int left_child;
    int right_child;
    int is_leaf;
    int prediction;
} TreeNode;

typedef struct {
    TreeNode* nodes;
    int count;
    int capacity;
} Tree;

typedef struct {
    float value;
    int label;
} FeatureLabelPair;

// --- Global Data ---
float* g_features;          // Dataset features: [num_samples][num_features]
int* g_labels;              // Dataset labels: [num_samples]
Tree* g_forest;             // Array of trees
volatile long long g_result; // Accumulated result for anti-optimization

// --- Workspace for computation (to avoid mallocs in timed section) ---
int* g_sample_indices_ws;   // For bootstrapping
int* g_feature_indices_ws;  // For random feature selection
FeatureLabelPair* g_pair_ws; // For sorting and finding best split
int* g_left_indices_ws;     // For partitioning
int* g_right_indices_ws;    // For partitioning

// --- Function Prototypes ---
void setup_benchmark(int argc, char *argv[]);
void run_computation();
void cleanup();
int build_tree_recursive(Tree* tree, int* sample_indices, int num_node_samples, int depth);
int make_leaf_node(Tree* tree, int* sample_indices, int num_node_samples);
int compare_pairs(const void* a, const void* b);

// --- Implementation ---
int compare_pairs(const void* a, const void* b) {
    float fa = ((FeatureLabelPair*)a)->value;
    float fb = ((FeatureLabelPair*)b)->value;
    if (fa < fb) return -1;
    if (fa > fb) return 1;
    return 0;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s num_samples num_features num_trees max_depth seed\n", argv[0]);
        exit(1);
    }
    P_NUM_SAMPLES = atoi(argv[1]);
    P_NUM_FEATURES = atoi(argv[2]);
    P_NUM_TREES = atoi(argv[3]);
    P_MAX_DEPTH = atoi(argv[4]);
    uint32_t seed = atoi(argv[5]);

    mt_seed(seed);

    P_FEATURES_PER_SPLIT = (int)sqrt((float)P_NUM_FEATURES);
    if (P_FEATURES_PER_SPLIT == 0) P_FEATURES_PER_SPLIT = 1;

    // Allocate data
    g_features = (float*)malloc((long)P_NUM_SAMPLES * P_NUM_FEATURES * sizeof(float));
    g_labels = (int*)malloc(P_NUM_SAMPLES * sizeof(int));
    g_forest = (Tree*)malloc(P_NUM_TREES * sizeof(Tree));
    if (!g_features || !g_labels || !g_forest) {
        fprintf(stderr, "Failed to allocate main data structures.\n");
        exit(1);
    }
    
    // Initialize data
    for (int i = 0; i < P_NUM_SAMPLES; i++) {
        for (int j = 0; j < P_NUM_FEATURES; j++) {
            g_features[i * P_NUM_FEATURES + j] = (float)mt_rand() / (float)UINT32_MAX;
        }
        g_labels[i] = mt_rand() % 2;
    }

    // Initialize trees with pre-allocated node arrays
    int max_nodes_per_tree = (1 << (P_MAX_DEPTH + 1)) - 1;
    for (int i = 0; i < P_NUM_TREES; i++) {
        g_forest[i].nodes = (TreeNode*)malloc(max_nodes_per_tree * sizeof(TreeNode));
        if (!g_forest[i].nodes) {
            fprintf(stderr, "Failed to allocate memory for tree nodes.\n");
            exit(1);
        }
        g_forest[i].count = 0;
        g_forest[i].capacity = max_nodes_per_tree;
    }

    // Allocate all workspace buffers
    g_sample_indices_ws = (int*)malloc(P_NUM_SAMPLES * sizeof(int));
    g_feature_indices_ws = (int*)malloc(P_NUM_FEATURES * sizeof(int));
    g_pair_ws = (FeatureLabelPair*)malloc(P_NUM_SAMPLES * sizeof(FeatureLabelPair));
    g_left_indices_ws = (int*)malloc(P_NUM_SAMPLES * sizeof(int));
    g_right_indices_ws = (int*)malloc(P_NUM_SAMPLES * sizeof(int));
    if (!g_sample_indices_ws || !g_feature_indices_ws || !g_pair_ws || !g_left_indices_ws || !g_right_indices_ws) {
        fprintf(stderr, "Failed to allocate workspace memory.\n");
        exit(1);
    }

    g_result = 0;
}

int make_leaf_node(Tree* tree, int* sample_indices, int num_node_samples) {
    if (tree->count >= tree->capacity) return -1; // Should not happen with pre-allocation
    int node_idx = tree->count++;
    TreeNode* node = &tree->nodes[node_idx];

    node->is_leaf = 1;
    node->left_child = -1;
    node->right_child = -1;
    node->feature_index = -1;
    node->threshold = 0.0f;

    if (num_node_samples == 0) {
        node->prediction = 0; // Default prediction
    } else {
        int count1 = 0;
        for (int i = 0; i < num_node_samples; i++) {
            if (g_labels[sample_indices[i]] == 1) {
                count1++;
            }
        }
        node->prediction = (count1 > num_node_samples / 2) ? 1 : 0;
    }
    return node_idx;
}

int build_tree_recursive(Tree* tree, int* sample_indices, int num_node_samples, int depth) {
    int first_label = -1;
    int is_pure = 1;
    if (num_node_samples > 0) {
        first_label = g_labels[sample_indices[0]];
        for (int i = 1; i < num_node_samples; i++) {
            if (g_labels[sample_indices[i]] != first_label) {
                is_pure = 0;
                break;
            }
        }
    }

    if (num_node_samples <= 1 || depth >= P_MAX_DEPTH || is_pure) {
        return make_leaf_node(tree, sample_indices, num_node_samples);
    }

    int best_feature = -1;
    float best_threshold = 0.0f;
    float best_gini = 1.1f; // Gini is always <= 1.0

    for (int i = 0; i < P_NUM_FEATURES; ++i) g_feature_indices_ws[i] = i;
    for (int i = P_NUM_FEATURES - 1; i > 0; --i) {
        int j = mt_rand() % (i + 1);
        int temp = g_feature_indices_ws[i]; g_feature_indices_ws[i] = g_feature_indices_ws[j]; g_feature_indices_ws[j] = temp;
    }
    
    int parent_counts[2] = {0, 0};
    for(int i = 0; i < num_node_samples; i++) parent_counts[g_labels[sample_indices[i]]]++;

    for (int f_idx = 0; f_idx < P_FEATURES_PER_SPLIT; f_idx++) {
        int feature = g_feature_indices_ws[f_idx];
        
        for (int i = 0; i < num_node_samples; i++) {
            g_pair_ws[i].value = g_features[sample_indices[i] * P_NUM_FEATURES + feature];
            g_pair_ws[i].label = g_labels[sample_indices[i]];
        }
        
        qsort(g_pair_ws, num_node_samples, sizeof(FeatureLabelPair), compare_pairs);
        
        int left_counts[2] = {0, 0};
        int right_counts[2] = {parent_counts[0], parent_counts[1]};

        for (int i = 0; i < num_node_samples - 1; i++) {
            int current_label = g_pair_ws[i].label;
            left_counts[current_label]++;
            right_counts[current_label]--;

            if (g_pair_ws[i].value == g_pair_ws[i+1].value) continue;

            float w_left = (float)(i + 1) / num_node_samples;
            float w_right = (float)(num_node_samples - i - 1) / num_node_samples;

            float gini_left = 1.0f - powf((float)left_counts[0] / (i + 1), 2) - powf((float)left_counts[1] / (i + 1), 2);
            float gini_right = (num_node_samples - i - 1 == 0) ? 0.0f : (1.0f - powf((float)right_counts[0] / (num_node_samples - i - 1), 2) - powf((float)right_counts[1] / (num_node_samples - i - 1), 2));
            
            float gini = w_left * gini_left + w_right * gini_right;

            if (gini < best_gini) {
                best_gini = gini;
                best_feature = feature;
                best_threshold = (g_pair_ws[i].value + g_pair_ws[i+1].value) / 2.0f;
            }
        }
    }
    
    if (best_feature == -1) {
        return make_leaf_node(tree, sample_indices, num_node_samples);
    }

    if (tree->count >= tree->capacity) return -1;
    int node_idx = tree->count++;
    TreeNode* node = &tree->nodes[node_idx];
    node->is_leaf = 0;
    node->feature_index = best_feature;
    node->threshold = best_threshold;
    node->prediction = -1;

    int left_count = 0, right_count = 0;
    for (int i = 0; i < num_node_samples; i++) {
        if (g_features[sample_indices[i] * P_NUM_FEATURES + best_feature] < best_threshold) {
             g_left_indices_ws[left_count++] = sample_indices[i];
        } else {
             g_right_indices_ws[right_count++] = sample_indices[i];
        }
    }
    
    g_result += node->feature_index;

    node->left_child = build_tree_recursive(tree, g_left_indices_ws, left_count, depth + 1);
    node->right_child = build_tree_recursive(tree, g_right_indices_ws, right_count, depth + 1);
    
    return node_idx;
}

void run_computation() {
    int* initial_indices = (int*)malloc(P_NUM_SAMPLES * sizeof(int));
    if (!initial_indices) { exit(1); } // Should have enough memory, but check
    for (int i = 0; i < P_NUM_SAMPLES; i++) {
        initial_indices[i] = i;
    }

    for (int i = 0; i < P_NUM_TREES; i++) {
        for (int j = 0; j < P_NUM_SAMPLES; j++) {
            g_sample_indices_ws[j] = initial_indices[mt_rand() % P_NUM_SAMPLES];
        }
        
        g_forest[i].count = 0;
        build_tree_recursive(&g_forest[i], g_sample_indices_ws, P_NUM_SAMPLES, 0);
        
        for(int k=0; k < g_forest[i].count; ++k) {
            if(g_forest[i].nodes[k].is_leaf){
                g_result += g_forest[i].nodes[k].prediction;
            }
        }
    }

    free(initial_indices);
}

void cleanup() {
    free(g_features);
    free(g_labels);
    
    for(int i = 0; i < P_NUM_TREES; i++) {
        free(g_forest[i].nodes);
    }
    free(g_forest);
    
    free(g_sample_indices_ws);
    free(g_feature_indices_ws);
    free(g_pair_ws);
    free(g_left_indices_ws);
    free(g_right_indices_ws);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
    cleanup();

    printf("%lld\n", g_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
