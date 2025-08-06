#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

// --- BEGIN MERSENNE TWISTER (DO NOT MODIFY) ---
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

// Benchmark parameters
static int num_samples;
static int num_features;
static int num_trees;
static int max_depth;
static float learning_rate;

// Data arrays
static float* features;       // Input features (X)
static float* targets;        // True values (y)
static float* predictions;    // Current model predictions
static float* residuals;      // Error = targets - predictions
static float* tree_predictions; // Predictions of the current tree being built

// Final result accumulator
static double final_result_accumulator = 0.0;

// Helper to generate a random float in [0, 1]
float rand_float() {
    return (float)mt_rand() / (float)UINT32_MAX;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 7) {
        fprintf(stderr, "Usage: %s num_samples num_features num_trees max_depth learning_rate seed\n", argv[0]);
        exit(1);
    }

    num_samples = atoi(argv[1]);
    num_features = atoi(argv[2]);
    num_trees = atoi(argv[3]);
    max_depth = atoi(argv[4]);
    learning_rate = atof(argv[5]);
    uint32_t seed = (uint32_t)atoi(argv[6]);

    mt_seed(seed);

    // Allocate memory
    features = (float*)malloc(num_samples * num_features * sizeof(float));
    targets = (float*)malloc(num_samples * sizeof(float));
    predictions = (float*)malloc(num_samples * sizeof(float));
    residuals = (float*)malloc(num_samples * sizeof(float));
    tree_predictions = (float*)malloc(num_samples * sizeof(float));

    if (!features || !targets || !predictions || !residuals || !tree_predictions) {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize data
    for (int i = 0; i < num_samples; ++i) {
        for (int j = 0; j < num_features; ++j) {
            features[i * num_features + j] = rand_float();
        }
        targets[i] = rand_float() * 10.0f; // Target values between 0 and 10
    }

    // Initial prediction is the mean of the targets
    double sum_targets = 0.0;
    for (int i = 0; i < num_samples; ++i) {
        sum_targets += targets[i];
    }
    float initial_prediction = (float)(sum_targets / num_samples);
    for (int i = 0; i < num_samples; ++i) {
        predictions[i] = initial_prediction;
    }
}

void run_computation() {
    // Main gradient boosting loop
    for (int t = 0; t < num_trees; ++t) {
        // Calculate residuals (negative gradient)
        for (int i = 0; i < num_samples; ++i) {
            residuals[i] = targets[i] - predictions[i];
        }

        // Reset predictions for the new tree
        memset(tree_predictions, 0, num_samples * sizeof(float));

        // Simulate training a regression tree on the residuals
        // This is a simplified process. At each level of depth, we find the single
        // best feature to split on among all samples and update leaf values.
        for (int d = 0; d < max_depth; ++d) {
            int best_feature = -1;
            float max_gain = -1.0f;

            // Find the best feature that reduces the error the most
            for (int f = 0; f < num_features; ++f) {
                float current_gain_left = 0.0f;
                float current_gain_right = 0.0f;
                 // A fixed split threshold is used for simplicity
                float threshold = 0.5f;

                for (int i = 0; i < num_samples; ++i) {
                    if (features[i * num_features + f] < threshold) {
                        current_gain_left += residuals[i];
                    } else {
                        current_gain_right += residuals[i];
                    }
                }
                // Simplified gain: absolute sum of residuals in partitions
                float total_gain = (current_gain_left > 0 ? current_gain_left : -current_gain_left) +
                                   (current_gain_right > 0 ? current_gain_right : -current_gain_right);

                if (total_gain > max_gain) {
                    max_gain = total_gain;
                    best_feature = f;
                }
            }

             // If a useful feature was found, update the tree's predictions
            if (best_feature != -1) {
                float threshold = 0.5f;
                float leaf_value_left = 0.0f, leaf_value_right = 0.0f;
                int count_left = 0, count_right = 0;

                for (int i = 0; i < num_samples; ++i) {
                    if (features[i * num_features + best_feature] < threshold) {
                        leaf_value_left += residuals[i];
                        count_left++;
                    } else {
                        leaf_value_right += residuals[i];
                        count_right++;
                    }
                }
                if (count_left > 0) leaf_value_left /= count_left;
                if (count_right > 0) leaf_value_right /= count_right;

                // Add this level's simple split prediction to the current tree
                for (int i = 0; i < num_samples; ++i) {
                    if (features[i * num_features + best_feature] < threshold) {
                        tree_predictions[i] += leaf_value_left;
                    } else {
                        tree_predictions[i] += leaf_value_right;
                    }
                }
            }
        }

        // Update the main model predictions with the new tree's output
        for (int i = 0; i < num_samples; ++i) {
            predictions[i] += learning_rate * tree_predictions[i];
        }
    }

    // Accumulate a final result to prevent dead code elimination.
    for (int i = 0; i < num_samples; ++i) {
        final_result_accumulator += predictions[i];
    }
}

void cleanup() {
    free(features);
    free(targets);
    free(predictions);
    free(residuals);
    free(tree_predictions);
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
    printf("%f\n", final_result_accumulator);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
