#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- Mersenne Twister (DO NOT MODIFY) ---
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
// --- End of Mersenne Twister ---

// --- Benchmark Globals ---
typedef struct {
    int num_samples;
    int num_features;
    int num_classes;
    float* features;  // Flattened 2D array: [sample][feature]
    int* labels;      // 1D array of class labels
    double final_result;
} BenchmarkData;

static BenchmarkData g_data;

// --- Helper Function ---

/**
 * @brief Calculates the Shannon entropy for a given subset of data.
 * Entropy H(S) = - sum_{i=1 to C} (p_i * log2(p_i))
 * @param sample_indices An array of indices for the samples in the subset.
 * @param count The number of samples in the subset.
 * @return The calculated entropy as a double.
 */
static double calculate_entropy(const int* sample_indices, int count) {
    if (count == 0) {
        return 0.0;
    }

    int* class_counts = (int*)calloc(g_data.num_classes, sizeof(int));
    if (!class_counts) {
        perror("Failed to allocate memory for class counts");
        exit(EXIT_FAILURE);
    }

    // Count occurrences of each class in the subset
    for (int i = 0; i < count; ++i) {
        int sample_idx = sample_indices[i];
        class_counts[g_data.labels[sample_idx]]++;
    }

    double entropy = 0.0;
    for (int i = 0; i < g_data.num_classes; ++i) {
        if (class_counts[i] > 0) {
            double p = (double)class_counts[i] / count;
            entropy -= p * log2(p);
        }
    }

    free(class_counts);
    return entropy;
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_samples> <num_features> <num_classes> <seed>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    g_data.num_samples = atoi(argv[1]);
    g_data.num_features = atoi(argv[2]);
    g_data.num_classes = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    if (g_data.num_samples <= 0 || g_data.num_features <= 0 || g_data.num_classes <= 1) {
        fprintf(stderr, "FATAL: Invalid parameters. Samples, features > 0 and classes > 1.\n");
        exit(EXIT_FAILURE);
    }
    
    mt_seed(seed);
    
    g_data.features = (float*)malloc((size_t)g_data.num_samples * g_data.num_features * sizeof(float));
    g_data.labels = (int*)malloc((size_t)g_data.num_samples * sizeof(int));
    g_data.final_result = 0.0;

    if (!g_data.features || !g_data.labels) {
        perror("Failed to allocate data memory");
        exit(EXIT_FAILURE);
    }

    // Generate random feature data (continuous values)
    for (int i = 0; i < g_data.num_samples; ++i) {
        for (int j = 0; j < g_data.num_features; ++j) {
            g_data.features[i * g_data.num_features + j] = (float)mt_rand() / (float)UINT32_MAX;
        }
    }

    // Generate random class labels (discrete values)
    for (int i = 0; i < g_data.num_samples; ++i) {
        g_data.labels[i] = mt_rand() % g_data.num_classes;
    }
}

void run_computation() {
    // Calculate the initial entropy of the entire dataset
    int* all_indices = (int*)malloc((size_t)g_data.num_samples * sizeof(int));
    if (!all_indices) {
        perror("Failed to allocate all_indices memory");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < g_data.num_samples; ++i) {
        all_indices[i] = i;
    }
    const double total_entropy = calculate_entropy(all_indices, g_data.num_samples);
    free(all_indices);

    // Pre-allocate memory for splitting indices to avoid malloc in the loop
    int* left_indices = (int*)malloc((size_t)g_data.num_samples * sizeof(int));
    int* right_indices = (int*)malloc((size_t)g_data.num_samples * sizeof(int));
    if (!left_indices || !right_indices) {
        perror("Failed to allocate split indices memory");
        exit(EXIT_FAILURE);
    }

    double total_gain_checksum = 0.0;

    // Iterate through each feature to find the best one to split on
    for (int f_idx = 0; f_idx < g_data.num_features; ++f_idx) {
        // Use the average value of the feature as the split threshold
        double threshold = 0.0;
        for (int s_idx = 0; s_idx < g_data.num_samples; ++s_idx) {
            threshold += g_data.features[s_idx * g_data.num_features + f_idx];
        }
        threshold /= g_data.num_samples;

        // Split the data into two subsets based on the threshold
        int left_count = 0;
        int right_count = 0;
        for (int s_idx = 0; s_idx < g_data.num_samples; ++s_idx) {
            if (g_data.features[s_idx * g_data.num_features + f_idx] <= threshold) {
                left_indices[left_count++] = s_idx;
            } else {
                right_indices[right_count++] = s_idx;
            }
        }

        // If a split results in an empty subset, it provides no information gain.
        if (left_count == 0 || right_count == 0) {
            continue;
        }
        
        // Calculate entropy of the two resulting subsets
        double left_entropy = calculate_entropy(left_indices, left_count);
        double right_entropy = calculate_entropy(right_indices, right_count);

        // Calculate the weighted average entropy of the split
        double p_left = (double)left_count / g_data.num_samples;
        double p_right = (double)right_count / g_data.num_samples;
        double weighted_entropy = (p_left * left_entropy) + (p_right * right_entropy);

        // Information gain is the reduction in entropy
        double gain = total_entropy - weighted_entropy;

        // Accumulate gain to prevent dead code elimination and have a meaningful result
        total_gain_checksum += gain;
    }
    
    g_data.final_result = total_gain_checksum;

    free(left_indices);
    free(right_indices);
}

void cleanup() {
    if (g_data.features) {
        free(g_data.features);
        g_data.features = NULL;
    }
    if (g_data.labels) {
        free(g_data.labels);
        g_data.labels = NULL;
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
    printf("%.6f\n", g_data.final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
