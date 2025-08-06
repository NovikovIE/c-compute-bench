#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- Mersenne Twister (MT19937) BEGIN ---
// Do not modify this section.
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
// --- Mersenne Twister (MT19937) END ---

// Benchmark parameters
int num_points;
int ambient_dimension;
float epsilon_radius;
int max_simplicial_dimension;

// Data structures
float* points;           // Flat array for point cloud: num_points * ambient_dimension
char* adjacency_matrix; // Flat array for graph: num_points * num_points
long long final_result;  // Result of the computation

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s num_points ambient_dimension epsilon_radius max_simplicial_dimension seed\n", argv[0]);
        exit(1);
    }

    num_points = atoi(argv[1]);
    ambient_dimension = atoi(argv[2]);
    epsilon_radius = atof(argv[3]);
    max_simplicial_dimension = atoi(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    mt_seed(seed);

    // Allocate memory for the point cloud
    points = (float*)malloc(num_points * ambient_dimension * sizeof(float));
    if (!points) {
        fprintf(stderr, "Failed to allocate memory for points\n");
        exit(1);
    }

    // Generate random points in a unit hypercube [0, 1]^d
    for (int i = 0; i < num_points * ambient_dimension; ++i) {
        points[i] = (float)mt_rand() / (float)UINT32_MAX;
    }

    // Allocate and build the adjacency matrix for the Vietoris-Rips graph
    adjacency_matrix = (char*)malloc((size_t)num_points * num_points * sizeof(char));
    if (!adjacency_matrix) {
        fprintf(stderr, "Failed to allocate memory for adjacency matrix\n");
        exit(1);
    }

    float epsilon_sq = epsilon_radius * epsilon_radius;

    for (int i = 0; i < num_points; ++i) {
        for (int j = i; j < num_points; ++j) {
            if (i == j) {
                adjacency_matrix[i * num_points + j] = 1;
                continue;
            }

            float dist_sq = 0.0f;
            for (int k = 0; k < ambient_dimension; ++k) {
                float diff = points[i * ambient_dimension + k] - points[j * ambient_dimension + k];
                dist_sq += diff * diff;
            }

            if (dist_sq <= epsilon_sq) {
                adjacency_matrix[i * num_points + j] = 1;
                adjacency_matrix[j * num_points + i] = 1;
            } else {
                adjacency_matrix[i * num_points + j] = 0;
                adjacency_matrix[j * num_points + i] = 0;
            }
        }
    }
}

void count_recursive_cliques(int* current_clique, int depth, int start_point_idx, const int target_depth) {
    if (depth == target_depth) {
        final_result++;
        return;
    }
    
    // Pruning: if we don't have enough remaining vertices to form the clique, stop.
    if (start_point_idx + (target_depth - depth) > num_points) {
        return;
    }

    for (int i = start_point_idx; i < num_points; ++i) {
        int is_fully_connected = 1;
        // Check connectivity with already selected vertices in the current clique
        for (int j = 0; j < depth; ++j) {
            if (adjacency_matrix[i * num_points + current_clique[j]] == 0) {
                is_fully_connected = 0;
                break;
            }
        }

        if (is_fully_connected) {
            current_clique[depth] = i;
            count_recursive_cliques(current_clique, depth + 1, i + 1, target_depth);
        }
    }
}

void run_computation() {
    final_result = 0;
    if (max_simplicial_dimension < 0) return;
    
    // A k-simplex is a (k+1)-clique in the graph.
    const int clique_size = max_simplicial_dimension + 1;

    if (clique_size == 1) { // 0-simplices are just the points
        final_result = num_points;
        return;
    }

    int* current_clique = (int*)malloc(clique_size * sizeof(int));
    if (!current_clique) {
        fprintf(stderr, "Failed to allocate memory for clique buffer\n");
        exit(1);
    }
    
    count_recursive_cliques(current_clique, 0, 0, clique_size);

    free(current_clique);
}

void cleanup() {
    free(points);
    free(adjacency_matrix);
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
