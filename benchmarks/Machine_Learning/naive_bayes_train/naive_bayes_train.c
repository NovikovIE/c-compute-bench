#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
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
// --- End of Mersenne Twister ---

// Benchmark parameters
long num_samples;
int num_features;
int num_classes;

// Assumed max value for a discrete feature (value will be in [0, MAX_FEATURE_VALUE-1])
#define MAX_FEATURE_VALUE 10

// Data structures for training
int* features; // Flattened 2D array: num_samples * num_features
int* labels;   // 1D array: num_samples

// Model parameters (computed during training)
long* class_counts;    // 1D array: num_classes
long* feature_counts;  // Flattened 3D array: num_classes * num_features * MAX_FEATURE_VALUE

// Result for verification to prevent dead code elimination
long final_result;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_samples> <num_features> <num_classes> <seed>\n", argv[0]);
        exit(1);
    }

    num_samples = atol(argv[1]);
    num_features = atoi(argv[2]);
    num_classes = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    mt_seed(seed);

    // Allocate memory
    features = (int*)malloc(num_samples * num_features * sizeof(int));
    labels = (int*)malloc(num_samples * sizeof(int));
    class_counts = (long*)malloc(num_classes * sizeof(long));
    long feature_counts_size = (long)num_classes * num_features * MAX_FEATURE_VALUE;
    feature_counts = (long*)malloc(feature_counts_size * sizeof(long));

    if (!features || !labels || !class_counts || !feature_counts) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Generate synthetic data
    for (long i = 0; i < num_samples; ++i) {
        labels[i] = mt_rand() % num_classes;
        for (int j = 0; j < num_features; ++j) {
            features[i * num_features + j] = mt_rand() % MAX_FEATURE_VALUE;
        }
    }
}

void run_computation() {
    // Initialize count arrays to zero
    memset(class_counts, 0, num_classes * sizeof(long));
    long feature_counts_size = (long)num_classes * num_features * MAX_FEATURE_VALUE;
    memset(feature_counts, 0, feature_counts_size * sizeof(long));

    // Training: count occurrences of each feature value for each class
    for (long i = 0; i < num_samples; ++i) {
        int current_label = labels[i];
        class_counts[current_label]++;
        for (int j = 0; j < num_features; ++j) {
            int feature_value = features[i * num_features + j];
            long index = ((long)current_label * num_features + j) * MAX_FEATURE_VALUE + feature_value;
            feature_counts[index]++;
        }
    }

    // Calculate a checksum to prevent dead code elimination
    final_result = 0;
    for (long i = 0; i < feature_counts_size; ++i) {
        final_result += feature_counts[i] * (i % 101); // Arbitrary calculation
    }
}

void cleanup() {
    free(features);
    free(labels);
    free(class_counts);
    free(feature_counts);
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
