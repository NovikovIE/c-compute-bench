#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <string.h>

// Mersenne Twister (verbatim)
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

// Benchmark-specific data structures
typedef struct {
    int NUM_POINTS;
    int NUM_FEATURES;
    int NUM_CLUSTERS;
    int MAX_ITERATIONS;

    float* points;
    float* centroids;
    int* assignments;
    
    // Auxiliary arrays for computation
    float* new_centroids;
    int* cluster_counts;

    long long final_result;
} BenchmarkData;

BenchmarkData data;

// Function to generate a random float between 0.0 and 1.0
float rand_float() {
    return (float)mt_rand() / (float)UINT32_MAX;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s num_points num_features num_clusters max_iterations seed\n", argv[0]);
        exit(1);
    }

    data.NUM_POINTS = atoi(argv[1]);
    data.NUM_FEATURES = atoi(argv[2]);
    data.NUM_CLUSTERS = atoi(argv[3]);
    data.MAX_ITERATIONS = atoi(argv[4]);
    uint32_t seed = (uint32_t)strtoul(argv[5], NULL, 10);
    
    mt_seed(seed);

    // Allocate memory
    data.points = (float*)malloc(data.NUM_POINTS * data.NUM_FEATURES * sizeof(float));
    data.centroids = (float*)malloc(data.NUM_CLUSTERS * data.NUM_FEATURES * sizeof(float));
    data.assignments = (int*)malloc(data.NUM_POINTS * sizeof(int));
    data.new_centroids = (float*)malloc(data.NUM_CLUSTERS * data.NUM_FEATURES * sizeof(float));
    data.cluster_counts = (int*)malloc(data.NUM_CLUSTERS * sizeof(int));
    
    if (!data.points || !data.centroids || !data.assignments || !data.new_centroids || !data.cluster_counts) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Generate random data points
    for (int i = 0; i < data.NUM_POINTS; ++i) {
        for (int j = 0; j < data.NUM_FEATURES; ++j) {
            data.points[i * data.NUM_FEATURES + j] = rand_float();
        }
    }

    // Initialize centroids using the Forgy method (select first K points)
    for (int i = 0; i < data.NUM_CLUSTERS; ++i) {
        for (int j = 0; j < data.NUM_FEATURES; ++j) {
            data.centroids[i * data.NUM_FEATURES + j] = data.points[i * data.NUM_FEATURES + j];
        }
    }
}

// Calculate squared Euclidean distance to avoid sqrt
float squared_distance(const float* p1, const float* p2, int features) {
    float dist = 0.0f;
    for (int i = 0; i < features; ++i) {
        float diff = p1[i] - p2[i];
        dist += diff * diff;
    }
    return dist;
}

void run_computation() {
    for (int iter = 0; iter < data.MAX_ITERATIONS; ++iter) {
        // --- Assignment Step ---
        // For each point, find the closest centroid
        for (int i = 0; i < data.NUM_POINTS; ++i) {
            float min_dist = FLT_MAX;
            int best_cluster = -1;
            const float* current_point = &data.points[i * data.NUM_FEATURES];

            for (int j = 0; j < data.NUM_CLUSTERS; ++j) {
                const float* current_centroid = &data.centroids[j * data.NUM_FEATURES];
                float dist = squared_distance(current_point, current_centroid, data.NUM_FEATURES);
                if (dist < min_dist) {
                    min_dist = dist;
                    best_cluster = j;
                }
            }
            data.assignments[i] = best_cluster;
        }

        // --- Update Step ---
        // Reset auxiliary arrays
        memset(data.new_centroids, 0, data.NUM_CLUSTERS * data.NUM_FEATURES * sizeof(float));
        memset(data.cluster_counts, 0, data.NUM_CLUSTERS * sizeof(int));

        // Sum up points for each cluster
        for (int i = 0; i < data.NUM_POINTS; ++i) {
            int cluster_idx = data.assignments[i];
            data.cluster_counts[cluster_idx]++;
            for (int j = 0; j < data.NUM_FEATURES; ++j) {
                data.new_centroids[cluster_idx * data.NUM_FEATURES + j] += data.points[i * data.NUM_FEATURES + j];
            }
        }
        
        // Calculate new centroids by averaging
        for (int i = 0; i < data.NUM_CLUSTERS; ++i) {
            if (data.cluster_counts[i] > 0) {
                for (int j = 0; j < data.NUM_FEATURES; ++j) {
                    data.centroids[i * data.NUM_FEATURES + j] = data.new_centroids[i * data.NUM_FEATURES + j] / data.cluster_counts[i];
                }
            }
            // else: cluster is empty, centroid remains unchanged from previous iteration.
        }
    }
    
    // Calculate final result to prevent dead code elimination
    long long total_assignments = 0;
    for(int i = 0; i < data.NUM_POINTS; ++i) {
        total_assignments += data.assignments[i];
    }
    data.final_result = total_assignments;
}

void cleanup() {
    free(data.points);
    free(data.centroids);
    free(data.assignments);
    free(data.new_centroids);
    free(data.cluster_counts);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    // Capture the result BEFORE cleanup to avoid use-after-free
    long long final_result = data.final_result;

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
    // Print result to stdout
    printf("%lld\n", final_result); 

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
