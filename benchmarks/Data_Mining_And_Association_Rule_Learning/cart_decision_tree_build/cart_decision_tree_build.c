#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- Mersenne Twister (MT19937) PRNG --- (DO NOT MODIFY)
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
// --- End of MT19937 --- 

// --- Benchmark Globals ---
int num_samples;
int num_features;
float *X_data; // Feature data (flattened 2D array)
int *y_data;   // Target labels (0 or 1)
double final_result; // Accumulated result to prevent dead code elimination

// --- CART-specific Helper Structures and Functions ---

// Struct to hold a feature value and its corresponding class label
typedef struct {
    float value;
    int label;
} ValueLabelPair;

// Comparison function for qsort to sort pairs by feature value
int compare_pairs(const void *a, const void *b) {
    ValueLabelPair *pair_a = (ValueLabelPair *)a;
    ValueLabelPair *pair_b = (ValueLabelPair *)b;
    if (pair_a->value < pair_b->value) return -1;
    if (pair_a->value > pair_b->value) return 1;
    return 0;
}

// Gini impurity calculation for a node/partition
// G = 1 - sum(p_i^2)
static inline float calculate_gini(int total_count, int class1_count) {
    if (total_count == 0) {
        return 0.0f;
    }
    float p1 = (float)class1_count / total_count;
    float p0 = 1.0f - p1;
    return 1.0f - (p0 * p0 + p1 * p1);
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_samples> <num_features> <seed>\n", argv[0]);
        exit(1);
    }

    num_samples = atoi(argv[1]);
    num_features = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (num_samples <= 0 || num_features <= 0) {
        fprintf(stderr, "FATAL: num_samples and num_features must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    X_data = (float *)malloc((size_t)num_samples * num_features * sizeof(float));
    y_data = (int *)malloc((size_t)num_samples * sizeof(int));

    if (!X_data || !y_data) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Generate random feature data (floats between 0.0 and 1.0)
    for (int i = 0; i < num_samples * num_features; ++i) {
        X_data[i] = (float)mt_rand() / (float)UINT32_MAX;
    }

    // Generate random binary class labels (0 or 1)
    for (int i = 0; i < num_samples; ++i) {
        y_data[i] = mt_rand() % 2;
    }

    final_result = 0.0;
}

void run_computation() {
    double total_min_gini = 0.0;

    // Allocate temporary storage for sorting feature values
    ValueLabelPair *feature_values = (ValueLabelPair *)malloc((size_t)num_samples * sizeof(ValueLabelPair));
    if (!feature_values) {
        fprintf(stderr, "FATAL: Temp memory allocation failed in computation.\n");
        exit(1);
    }

    // Calculate total number of samples in class 1 initially
    int total_class1_samples = 0;
    for (int i = 0; i < num_samples; ++i) {
        if (y_data[i] == 1) {
            total_class1_samples++;
        }
    }

    // Iterate over each feature to find the best split point
    for (int f = 0; f < num_features; ++f) {
        // Populate the temporary array with values and labels for the current feature
        for (int i = 0; i < num_samples; ++i) {
            feature_values[i].value = X_data[i * num_features + f];
            feature_values[i].label = y_data[i];
        }

        // Sort the samples based on the current feature's value
        qsort(feature_values, num_samples, sizeof(ValueLabelPair), compare_pairs);

        float min_gini_for_feature = 1.0f;

        // Initialize counts for left and right partitions
        int left_count = 0;
        int left_class1_count = 0;
        int right_count = num_samples;
        int right_class1_count = total_class1_samples;

        // Iterate through all possible split points (between sorted values)
        for (int i = 0; i < num_samples - 1; ++i) {
            // Move one sample from the 'right' partition to the 'left'
            left_count++;
            right_count--;
            if (feature_values[i].label == 1) {
                left_class1_count++;
                right_class1_count--;
            }
            
            // If the next value is the same, don't evaluate a split here
            if (feature_values[i].value == feature_values[i+1].value) {
                continue;
            }

            // Calculate Gini impurity for both partitions
            float gini_left = calculate_gini(left_count, left_class1_count);
            float gini_right = calculate_gini(right_count, right_class1_count);

            // Calculate weighted average of Gini impurities
            float weighted_gini = ((float)left_count / num_samples) * gini_left + 
                                  ((float)right_count / num_samples) * gini_right;

            if (weighted_gini < min_gini_for_feature) {
                min_gini_for_feature = weighted_gini;
            }
        }

        // Accumulate the best Gini found for this feature
        total_min_gini += min_gini_for_feature;
    }

    free(feature_values);
    final_result = total_min_gini;
}

void cleanup() {
    free(X_data);
    free(y_data);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final accumulated result to stdout
    printf("%f\n", final_result);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
