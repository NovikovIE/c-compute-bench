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

// --- Benchmark Globals ---
int num_points;
int num_features;
float** points;
float** dist_matrix;
int* cluster_map;
float final_result; // Accumulator for preventing dead code elimination

// --- Benchmark Functions ---

/**
 * @brief Sets up the benchmark data.
 * Parses command-line arguments for num_points, num_features, and seed.
 * Allocates memory for data points, a distance matrix, and cluster mappings.
 * Generates random data points and pre-computes the pairwise Euclidean distance matrix.
 */
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s num_points num_features seed\n", argv[0]);
        exit(1);
    }

    num_points = atoi(argv[1]);
    num_features = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (num_points <= 1) {
        fprintf(stderr, "Error: num_points must be greater than 1.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate memory for data points
    points = (float**)malloc(num_points * sizeof(float*));
    for (int i = 0; i < num_points; i++) {
        points[i] = (float*)malloc(num_features * sizeof(float));
    }

    // Allocate memory for the distance matrix and cluster map
    dist_matrix = (float**)malloc(num_points * sizeof(float*));
    for (int i = 0; i < num_points; i++) {
        dist_matrix[i] = (float*)malloc(num_points * sizeof(float));
    }
    cluster_map = (int*)malloc(num_points * sizeof(int));

    // Generate random data points in the range [0.0, 1.0]
    for (int i = 0; i < num_points; i++) {
        for (int j = 0; j < num_features; j++) {
            points[i][j] = (float)mt_rand() / (float)UINT32_MAX;
        }
    }

    // Pre-compute the pairwise Euclidean distance matrix
    for (int i = 0; i < num_points; i++) {
        dist_matrix[i][i] = FLT_MAX; // Distance to self is effectively infinite
        for (int j = i + 1; j < num_points; j++) {
            float sum_sq_diff = 0.0f;
            for (int k = 0; k < num_features; k++) {
                float diff = points[i][k] - points[j][k];
                sum_sq_diff += diff * diff;
            }
            float dist = sqrtf(sum_sq_diff);
            dist_matrix[i][j] = dist;
            dist_matrix[j][i] = dist; // Matrix is symmetric
        }
    }
}

/**
 * @brief Runs the core computation of the benchmark.
 * Performs a naive agglomerative hierarchical clustering.
 * It iteratively finds the two closest clusters and merges them.
 * The sum of the distances at each merge is accumulated into `final_result`.
 * This is a O(N^3) implementation, suitable for a CPU-intensive benchmark.
 */
void run_computation() {
    float total_distance_sum = 0.0f;

    // Initialize: each point is its own cluster representative
    for (int i = 0; i < num_points; i++) {
        cluster_map[i] = i;
    }

    // Perform (num_points - 1) merges to form a single cluster
    for (int k = 0; k < num_points - 1; k++) {
        float min_dist = FLT_MAX;
        int merge_c1 = -1;
        int merge_c2 = -1;

        // Find the two closest active clusters
        // An active cluster is one which is its own representative (cluster_map[i] == i)
        for (int i = 0; i < num_points; i++) {
            if (cluster_map[i] != i) continue; // Skip if not a cluster representative

            for (int j = i + 1; j < num_points; j++) {
                if (cluster_map[j] != j) continue; // Skip if not a cluster representative

                if (dist_matrix[i][j] < min_dist) {
                    min_dist = dist_matrix[i][j];
                    merge_c1 = i;
                    merge_c2 = j;
                }
            }
        }

        if (merge_c1 != -1 && merge_c2 != -1) {
            // Merge the smaller index cluster into the larger one
            cluster_map[merge_c2] = merge_c1;
            total_distance_sum += min_dist;
        } else {
            // Should not happen in a valid run
            break;
        }
    }

    final_result = total_distance_sum;
}

/**
 * @brief Frees all memory allocated in setup_benchmark.
 */
void cleanup() {
    for (int i = 0; i < num_points; i++) {
        free(points[i]);
    }
    free(points);

    for (int i = 0; i < num_points; i++) {
        free(dist_matrix[i]);
    }
    free(dist_matrix);

    free(cluster_map);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final accumulated result to stdout
    printf("%f\n", final_result);

    // Print timing information to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
