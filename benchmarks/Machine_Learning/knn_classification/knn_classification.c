#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <float.h>

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
// --- End of Mersenne Twister ---

// --- Benchmark Globals ---
// Parameters
int num_training_points;
int num_test_points;
int num_features;
int k_neighbors;
const int num_classes = 10; // Fixed number of classes for simplicity

// Data
float* training_points;
int*   training_labels;
float* test_points;
int*   predictions;

// Result accumulator
long long final_result;

// Helper struct for storing neighbor information
typedef struct {
    float distance;
    int label;
} Neighbor;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <num_training_points> <num_test_points> <num_features> <k_neighbors> <seed>\n", argv[0]);
        exit(1);
    }

    num_training_points = atoi(argv[1]);
    num_test_points     = atoi(argv[2]);
    num_features        = atoi(argv[3]);
    k_neighbors         = atoi(argv[4]);
    uint32_t seed       = atoi(argv[5]);

    if (k_neighbors <= 0 || k_neighbors > num_training_points || num_classes <= 0) {
        fprintf(stderr, "FATAL: Invalid parameters. k > 0, k <= n_train, n_classes > 0\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate memory
    training_points = (float*)malloc((size_t)num_training_points * num_features * sizeof(float));
    training_labels = (int*)malloc((size_t)num_training_points * sizeof(int));
    test_points     = (float*)malloc((size_t)num_test_points * num_features * sizeof(float));
    predictions     = (int*)malloc((size_t)num_test_points * sizeof(int));

    if (!training_points || !training_labels || !test_points || !predictions) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize training data with random values
    for (int i = 0; i < num_training_points; ++i) {
        for (int j = 0; j < num_features; ++j) {
            training_points[i * num_features + j] = (float)mt_rand() / (float)UINT32_MAX;
        }
        training_labels[i] = mt_rand() % num_classes;
    }

    // Initialize test data with random values
    for (int i = 0; i < num_test_points; ++i) {
        for (int j = 0; j < num_features; ++j) {
            test_points[i * num_features + j] = (float)mt_rand() / (float)UINT32_MAX;
        }
    }
}

void run_computation() {
    // Allocate temporary arrays for neighbor tracking and voting
    Neighbor* k_nearest_neighbors = (Neighbor*)malloc(k_neighbors * sizeof(Neighbor));
    int* class_votes = (int*)malloc(num_classes * sizeof(int));
    if (!k_nearest_neighbors || !class_votes) {
        fprintf(stderr, "FATAL: Memory allocation failed for temporary arrays in computation.\n");
        exit(1);
    }

    // For each test point
    for (int i = 0; i < num_test_points; ++i) {
        // Initialize k-nearest neighbors with maximum distance
        for (int k = 0; k < k_neighbors; ++k) {
            k_nearest_neighbors[k].distance = FLT_MAX;
        }

        // For each training point, find its distance to the current test point
        for (int j = 0; j < num_training_points; ++j) {
            // Calculate squared Euclidean distance (avoids sqrt for comparison)
            float dist_sq = 0.0f;
            for (int f = 0; f < num_features; ++f) {
                float diff = test_points[i * num_features + f] - training_points[j * num_features + f];
                dist_sq += diff * diff;
            }

            // Find the neighbor with the largest distance in the current k-list
            int max_dist_idx = 0;
            for (int k = 1; k < k_neighbors; ++k) {
                if (k_nearest_neighbors[k].distance > k_nearest_neighbors[max_dist_idx].distance) {
                    max_dist_idx = k;
                }
            }

            // If the current point is closer, replace the farthest neighbor
            if (dist_sq < k_nearest_neighbors[max_dist_idx].distance) {
                k_nearest_neighbors[max_dist_idx].distance = dist_sq;
                k_nearest_neighbors[max_dist_idx].label = training_labels[j];
            }
        }

        // Perform majority vote among the k-nearest neighbors
        for(int c = 0; c < num_classes; ++c) {
            class_votes[c] = 0;
        }
        for (int k = 0; k < k_neighbors; ++k) {
            class_votes[k_nearest_neighbors[k].label]++;
        }

        int majority_class = 0;
        int max_votes = -1;
        for (int c = 0; c < num_classes; ++c) {
            if (class_votes[c] > max_votes) {
                max_votes = class_votes[c];
                majority_class = c;
            }
        }
        predictions[i] = majority_class;
    }

    // Accumulate result to prevent dead code elimination
    final_result = 0;
    for (int i = 0; i < num_test_points; ++i) {
        final_result += predictions[i];
    }
    
    free(k_nearest_neighbors);
    free(class_votes);
}

void cleanup() {
    free(training_points);
    free(training_labels);
    free(test_points);
    free(predictions);
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
    printf("%lld\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
