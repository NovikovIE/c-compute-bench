#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

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
// --- End of MT19937 ---

// Benchmark parameters and data structures
static int NUM_POINTS;
static int NUM_FEATURES;
static double EPSILON;
static int MIN_SAMPLES;

static double** points;
static int* cluster_ids; // 0: unclassified, -1: noise, >0: cluster ID

static int final_result;

// Helper to calculate squared Euclidean distance to avoid sqrt
static double euclidean_distance_sq(int p1_idx, int p2_idx) {
    double sum = 0.0;
    for (int i = 0; i < NUM_FEATURES; i++) {
        double diff = points[p1_idx][i] - points[p2_idx][i];
        sum += diff * diff;
    }
    return sum;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s num_points num_features epsilon min_samples seed\n", argv[0]);
        exit(1);
    }

    NUM_POINTS = atoi(argv[1]);
    NUM_FEATURES = atoi(argv[2]);
    EPSILON = atof(argv[3]);
    MIN_SAMPLES = atoi(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    mt_seed(seed);

    points = (double**)malloc(NUM_POINTS * sizeof(double*));
    cluster_ids = (int*)malloc(NUM_POINTS * sizeof(int));
    if (!points || !cluster_ids) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    for (int i = 0; i < NUM_POINTS; i++) {
        points[i] = (double*)malloc(NUM_FEATURES * sizeof(double));
        if (!points[i]) {
            fprintf(stderr, "Memory allocation failed for point %d\n", i);
            exit(1);
        }
        for (int j = 0; j < NUM_FEATURES; j++) {
            points[i][j] = ((double)mt_rand() / (double)UINT32_MAX) * 100.0;
        }
        cluster_ids[i] = 0; // 0 for unclassified
    }
}

void run_computation() {
    int cluster_id_counter = 0;
    int* queue = (int*)malloc(NUM_POINTS * sizeof(int));
    int* neighbors_buffer = (int*)malloc(NUM_POINTS * sizeof(int));
    if (!queue || !neighbors_buffer) {
        fprintf(stderr, "Buffer allocation failed in computation\n");
        exit(1);
    }

    double epsilon_sq = EPSILON * EPSILON;

    for (int i = 0; i < NUM_POINTS; i++) {
        if (cluster_ids[i] != 0) continue; // Already processed

        // Find neighbors of point i (Region Query)
        int num_neighbors = 0;
        for (int j = 0; j < NUM_POINTS; j++) {
            if (euclidean_distance_sq(i, j) <= epsilon_sq) {
                neighbors_buffer[num_neighbors++] = j;
            }
        }

        if (num_neighbors < MIN_SAMPLES) {
            cluster_ids[i] = -1; // Mark as noise
            continue;
        }

        // Core point found, start a new cluster
        cluster_id_counter++;
        int queue_head = 0;
        int queue_tail = 0;
        for (int k = 0; k < num_neighbors; k++) {
            queue[queue_tail++] = neighbors_buffer[k];
        }

        // Process the queue to expand the cluster
        while (queue_head < queue_tail) {
            int current_p_idx = queue[queue_head++];

            // Skip if already in a different cluster (can happen with border points of dense clusters)
            if (cluster_ids[current_p_idx] > 0 && cluster_ids[current_p_idx] != cluster_id_counter) continue;
            
            cluster_ids[current_p_idx] = cluster_id_counter;

            // If this point is also a core point, add its neighbors to the queue
            int inner_num_neighbors = 0;
            for (int j = 0; j < NUM_POINTS; j++) {
                if (euclidean_distance_sq(current_p_idx, j) <= epsilon_sq) {
                    inner_num_neighbors++;
                }
            }

            if (inner_num_neighbors >= MIN_SAMPLES) {
                for (int j = 0; j < NUM_POINTS; j++) {
                    if (euclidean_distance_sq(current_p_idx, j) <= epsilon_sq) {
                        if (cluster_ids[j] == 0) { // Add unclassified points to queue
                            cluster_ids[j] = cluster_id_counter;
                            if (queue_tail < NUM_POINTS) queue[queue_tail++] = j;
                        } else if (cluster_ids[j] == -1) { // Claim noise points
                            cluster_ids[j] = cluster_id_counter;
                        }
                    }
                }
            }
        }
    }

    free(queue);
    free(neighbors_buffer);

    // Use the result to prevent dead code elimination.
    // A simple checksum of cluster assignments.
    long long checksum = 0;
    for (int i = 0; i < NUM_POINTS; i++) {
        checksum += cluster_ids[i];
    }
    final_result = cluster_id_counter + (int)(checksum & 0xFFFFFF);
}

void cleanup() {
    for (int i = 0; i < NUM_POINTS; i++) {
        free(points[i]);
    }
    free(points);
    free(cluster_ids);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%d\n", final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
