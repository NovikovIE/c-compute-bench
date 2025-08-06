/**
 * adaboost_train.c
 * 
 * A benchmark simulating the training phase of an AdaBoost (Adaptive Boosting) classifier.
 * AdaBoost is an ensemble learning algorithm that combines multiple weak learners (in this case, decision stumps)
 * to create a strong classifier. This benchmark focuses on the iterative process of fitting weak learners
 * on re-weighted training data.
 *
 * The core computation involves iterating through a number of estimators (weak learners).
 * For each estimator, it finds the best decision stump by selecting a feature and a threshold
 * that best splits the data according to the current sample weights. A decision stump is a simple
 * model that classifies samples based on whether a single feature value is above or below a threshold.
 * After finding the best stump, the algorithm calculates its influence (alpha) and updates the sample weights,
 * giving more weight to previously misclassified samples. This forces the next estimator to focus
 * on the harder-to-classify examples.
 *
 * Parameters:
 * - num_samples: The number of data points in the training set.
 * - num_features: The number of features for each data point.
 * - num_estimators: The number of weak learners (decision stumps) to train and combine.
 * - seed: The seed for the random number generator.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>


// --- BEGIN MERSENNE TWISTER --- (Do Not Modify)
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

// Decision stump (weak learner) structure
typedef struct {
    int feature_index;  // Index of the feature to split on
    double threshold;   // Threshold for the feature
    double alpha;       // Weight of this stump in the final classifier
    int polarity;       // 1 or -1, to handle 'less than' or 'greater than' splits
} DecisionStump;

// Benchmark parameters
static int num_samples;
static int num_features;
static int num_estimators;

// Data arrays
static double** X_train;  // Feature data: num_samples x num_features
static int* y_train;      // Labels: num_samples

// AdaBoost-specific arrays
static double* sample_weights;
static DecisionStump* estimators;

// Final result accumulator
static double final_result = 0.0;

// Helper to generate a random double between 0.0 and 1.0
double rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_samples num_features num_estimators seed\n", argv[0]);
        exit(1);
    }

    num_samples = atoi(argv[1]);
    num_features = atoi(argv[2]);
    num_estimators = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    mt_seed(seed);

    // Allocate feature matrix X_train contiguously for better cache performance
    double* X_data = (double*)malloc(num_samples * num_features * sizeof(double));
    X_train = (double**)malloc(num_samples * sizeof(double*));
    if (!X_data || !X_train) {
        fprintf(stderr, "Memory allocation failed for X_train.\n"); exit(1);
    }
    for (int i = 0; i < num_samples; i++) {
        X_train[i] = &X_data[i * num_features];
    }

    // Allocate labels, weights, and estimators
    y_train = (int*)malloc(num_samples * sizeof(int));
    sample_weights = (double*)malloc(num_samples * sizeof(double));
    estimators = (DecisionStump*)malloc(num_estimators * sizeof(DecisionStump));
    if (!y_train || !sample_weights || !estimators) {
        fprintf(stderr, "Memory allocation failed for labels/weights/estimators.\n"); exit(1);
    }

    // Generate random training data
    for (int i = 0; i < num_samples; i++) {
        for (int j = 0; j < num_features; j++) {
            X_train[i][j] = rand_double();
        }
        // Generate binary labels: -1 or 1
        y_train[i] = (mt_rand() % 2) * 2 - 1;
    }
}

void run_computation() {
    // 1. Initialize weights for all samples equally
    for (int i = 0; i < num_samples; i++) {
        sample_weights[i] = 1.0 / num_samples;
    }

    final_result = 0.0;

    // 2. Main AdaBoost loop: iterate for each estimator
    for (int m = 0; m < num_estimators; m++) {
        double min_error = 1.0 / 0.0; // Infinity
        DecisionStump best_stump = {0};

        // Find the best decision stump (weak learner)
        for (int j = 0; j < num_features; j++) {
            // For this benchmark, we use the simple mean of the feature as a threshold
            double feature_sum = 0.0;
            for (int i = 0; i < num_samples; i++) {
                feature_sum += X_train[i][j];
            }
            double threshold = feature_sum / num_samples;

            for (int p = 0; p < 2; p++) {
                int polarity = (p == 0) ? 1 : -1;
                double error = 0.0;

                // Calculate weighted error of this stump
                for (int i = 0; i < num_samples; i++) {
                    int prediction = (X_train[i][j] < threshold) ? -1 : 1;
                    if (prediction * polarity != y_train[i]) {
                        error += sample_weights[i];
                    }
                }

                if (error < min_error) {
                    min_error = error;
                    best_stump.feature_index = j;
                    best_stump.threshold = threshold;
                    best_stump.polarity = polarity;
                }
            }
        }

        // 3. Calculate estimator weight (alpha)
        // Clamp error to avoid division by zero or log of non-positive
        double epsilon = fmax(min_error, 1e-10);
        best_stump.alpha = 0.5 * log((1.0 - epsilon) / epsilon);

        // 4. Update sample weights
        double weight_sum = 0.0;
        for (int i = 0; i < num_samples; i++) {
            int prediction = (X_train[i][best_stump.feature_index] < best_stump.threshold) ? -1 : 1;
            sample_weights[i] *= exp(-best_stump.alpha * y_train[i] * (prediction * best_stump.polarity));
            weight_sum += sample_weights[i];
        }

        // 5. Normalize weights
        for (int i = 0; i < num_samples; i++) {
            sample_weights[i] /= weight_sum;
        }

        // Store the trained stump
        estimators[m] = best_stump;
        final_result += best_stump.alpha; // Accumulate result
    }
}

void cleanup() {
    if (X_train) free(X_train[0]);
    free(X_train);
    free(y_train);
    free(sample_weights);
    free(estimators);
}

int main(int argc, char* argv[]) {
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
