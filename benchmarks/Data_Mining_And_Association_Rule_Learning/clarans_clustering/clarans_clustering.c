#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <string.h>

// --- MERSENNE TWISTER (DO NOT MODIFY) ---
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

// --- BENCHMARK DATA & PARAMETERS ---
typedef struct {
    int num_points;
    int num_features;
    int num_clusters;
    int num_local_searches;
    int max_neighbors;
    int num_random_numbers;

    float **points;
    int *current_medoids; 
    int *best_medoids;
    int *neighbor_medoids;
    uint32_t *random_numbers;

    double final_result;
    int rand_idx;
} BenchmarkData;

BenchmarkData g_data;

// --- FORWARD DECLARATIONS ---
float euclidean_dist_sq(int p_idx1, int p_idx2);
double calculate_total_cost(const int* medoids);

// --- DATA SETUP ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 8) {
        fprintf(stderr, "Usage: %s num_points num_features num_clusters num_local max_neighbors num_random seed\n", argv[0]);
        exit(1);
    }

    g_data.num_points = atoi(argv[1]);
    g_data.num_features = atoi(argv[2]);
    g_data.num_clusters = atoi(argv[3]);
    g_data.num_local_searches = atoi(argv[4]);
    g_data.max_neighbors = atoi(argv[5]);
    g_data.num_random_numbers = atoi(argv[6]);
    uint32_t seed = (uint32_t)atoi(argv[7]);
    
    mt_seed(seed);

    // Allocate data points
    g_data.points = (float **)malloc(g_data.num_points * sizeof(float *));
    for (int i = 0; i < g_data.num_points; ++i) {
        g_data.points[i] = (float *)malloc(g_data.num_features * sizeof(float));
        for (int j = 0; j < g_data.num_features; ++j) {
            g_data.points[i][j] = (float)mt_rand() / (float)UINT32_MAX;
        }
    }

    // Allocate medoid arrays
    g_data.current_medoids = (int *)malloc(g_data.num_clusters * sizeof(int));
    g_data.best_medoids = (int *)malloc(g_data.num_clusters * sizeof(int));
    g_data.neighbor_medoids = (int *)malloc(g_data.num_clusters * sizeof(int));
    
    // Pre-generate random numbers
    g_data.random_numbers = (uint32_t *) malloc(g_data.num_random_numbers * sizeof(uint32_t));
    for(int i = 0; i < g_data.num_random_numbers; ++i) {
        g_data.random_numbers[i] = mt_rand();
    }

    g_data.final_result = 0.0;
    g_data.rand_idx = 0;
}

// --- COMPUTATION ---
void run_computation() {
    double min_cost = INFINITY;
    g_data.rand_idx = 0;

    for (int i = 0; i < g_data.num_local_searches; ++i) {
        // 1. Select a random initial set of medoids for 'current_medoids'
        for (int k = 0; k < g_data.num_clusters; ++k) {
            int candidate;
            int is_duplicate;
            do {
                is_duplicate = 0;
                candidate = g_data.random_numbers[g_data.rand_idx++ % g_data.num_random_numbers] % g_data.num_points;
                for (int j = 0; j < k; ++j) {
                    if (g_data.current_medoids[j] == candidate) {
                        is_duplicate = 1;
                        break;
                    }
                }
            } while (is_duplicate);
            g_data.current_medoids[k] = candidate;
        }

        // 2. Calculate the cost of the current set
        double current_cost = calculate_total_cost(g_data.current_medoids);

        // 3. Explore neighbors
        for (int j = 0; j < g_data.max_neighbors; ++j) {
            // Create a neighbor by swapping one medoid with one non-medoid
            memcpy(g_data.neighbor_medoids, g_data.current_medoids, g_data.num_clusters * sizeof(int));

            int medoid_to_swap_idx = g_data.random_numbers[g_data.rand_idx++ % g_data.num_random_numbers] % g_data.num_clusters;

            int replacement_point_idx;
            int is_medoid;
            do {
                is_medoid = 0;
                replacement_point_idx = g_data.random_numbers[g_data.rand_idx++ % g_data.num_random_numbers] % g_data.num_points;
                for (int k = 0; k < g_data.num_clusters; ++k) {
                    if (g_data.current_medoids[k] == replacement_point_idx) {
                        is_medoid = 1;
                        break;
                    }
                }
            } while (is_medoid);

            g_data.neighbor_medoids[medoid_to_swap_idx] = replacement_point_idx;

            double neighbor_cost = calculate_total_cost(g_data.neighbor_medoids);

            if (neighbor_cost < current_cost) {
                current_cost = neighbor_cost;
                memcpy(g_data.current_medoids, g_data.neighbor_medoids, g_data.num_clusters * sizeof(int));
                // In the classic CLARANS, one might reset the neighbor search (j=0).
                // For this benchmark, we'll just continue the local search from the new, better point.
            }
        }

        if (current_cost < min_cost) {
            min_cost = current_cost;
            memcpy(g_data.best_medoids, g_data.current_medoids, g_data.num_clusters * sizeof(int));
        }
    }
    g_data.final_result = min_cost;
}


// --- HELPER FUNCTIONS ---
float euclidean_dist_sq(int p_idx1, int p_idx2) {
    float dist = 0.0f;
    float diff;
    for (int i = 0; i < g_data.num_features; ++i) {
        diff = g_data.points[p_idx1][i] - g_data.points[p_idx2][i];
        dist += diff * diff;
    }
    return dist;
}

double calculate_total_cost(const int* medoids) {
    double total_cost = 0.0;
    for (int i = 0; i < g_data.num_points; ++i) {
        float min_dist_sq = INFINITY;
        for (int j = 0; j < g_data.num_clusters; ++j) {
            float dist_sq = euclidean_dist_sq(i, medoids[j]);
            if (dist_sq < min_dist_sq) {
                min_dist_sq = dist_sq;
            }
        }
        total_cost += sqrt(min_dist_sq); // Use actual distance for total cost as is standard
    }
    return total_cost;
}

// --- CLEANUP ---
void cleanup() {
    for (int i = 0; i < g_data.num_points; ++i) {
        free(g_data.points[i]);
    }
    free(g_data.points);
    free(g_data.current_medoids);
    free(g_data.best_medoids);
    free(g_data.neighbor_medoids);
    free(g_data.random_numbers);
}

// --- MAIN ---
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%f\n", g_data.final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
