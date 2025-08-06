#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>

// --- BEGIN MERSENNE TWISTER (MT19937) --- (DO NOT MODIFY)
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

// Global struct to hold all benchmark data
typedef struct {
    int num_landmark_points;
    int num_witness_points;
    int max_simplicial_dimension;
    int point_dimension;
    double max_distance;
    uint32_t seed;

    double** landmarks;
    double** witnesses;

    unsigned long long simplex_count;
    double max_dist_sq; // Pre-calculated for efficiency

} BenchmarkData;

BenchmarkData g_data;

// Forward declarations for helper functions
bool has_witness(const int* simplex_indices, int simplex_size);
void count_simplices_recursive(int* simplex_indices, int current_size, int start_landmark_idx);

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 7) {
        fprintf(stderr, "Usage: %s num_landmark_points num_witness_points max_simplicial_dimension point_dimension max_distance seed\n", argv[0]);
        exit(1);
    }

    // Parse arguments
    g_data.num_landmark_points = atoi(argv[1]);
    g_data.num_witness_points = atoi(argv[2]);
    g_data.max_simplicial_dimension = atoi(argv[3]);
    g_data.point_dimension = atoi(argv[4]);
    g_data.max_distance = atof(argv[5]);
    g_data.seed = (uint32_t)strtoul(argv[6], NULL, 10);

    // Seed the random number generator
    mt_seed(g_data.seed);

    // Pre-calculate squared distance to avoid repeated sqrt in computation
    g_data.max_dist_sq = g_data.max_distance * g_data.max_distance;
    g_data.simplex_count = 0;

    // Allocate memory for landmark points
    g_data.landmarks = (double**)malloc(g_data.num_landmark_points * sizeof(double*));
    if (!g_data.landmarks) { perror("malloc failed"); exit(1); }
    for (int i = 0; i < g_data.num_landmark_points; ++i) {
        g_data.landmarks[i] = (double*)malloc(g_data.point_dimension * sizeof(double));
        if (!g_data.landmarks[i]) { perror("malloc failed"); exit(1); }
    }

    // Allocate memory for witness points
    g_data.witnesses = (double**)malloc(g_data.num_witness_points * sizeof(double*));
    if (!g_data.witnesses) { perror("malloc failed"); exit(1); }
    for (int i = 0; i < g_data.num_witness_points; ++i) {
        g_data.witnesses[i] = (double*)malloc(g_data.point_dimension * sizeof(double));
        if (!g_data.witnesses[i]) { perror("malloc failed"); exit(1); }
    }

    // Generate random landmark points in a unit hypercube [0, 1]^d
    for (int i = 0; i < g_data.num_landmark_points; ++i) {
        for (int j = 0; j < g_data.point_dimension; ++j) {
            g_data.landmarks[i][j] = (double)mt_rand() / (double)UINT32_MAX;
        }
    }

    // Generate random witness points
    for (int i = 0; i < g_data.num_witness_points; ++i) {
        for (int j = 0; j < g_data.point_dimension; ++j) {
            g_data.witnesses[i][j] = (double)mt_rand() / (double)UINT32_MAX;
        }
    }
}

// Checks if a given simplex (defined by a set of landmark indices) has a witness.
bool has_witness(const int* simplex_indices, int simplex_size) {
    for (int w_idx = 0; w_idx < g_data.num_witness_points; ++w_idx) {
        double max_dist_sq_for_this_witness = 0.0;
        // Find the max squared distance from this witness to any landmark in the simplex
        for (int s_idx = 0; s_idx < simplex_size; ++s_idx) {
            int l_idx = simplex_indices[s_idx];
            double dist_sq = 0.0;
            // Calculate squared Euclidean distance
            for (int d = 0; d < g_data.point_dimension; ++d) {
                double diff = g_data.witnesses[w_idx][d] - g_data.landmarks[l_idx][d];
                dist_sq += diff * diff;
            }
            if (dist_sq > max_dist_sq_for_this_witness) {
                max_dist_sq_for_this_witness = dist_sq;
            }
        }

        // If the max distance is within the threshold, we found a witness for this simplex
        if (max_dist_sq_for_this_witness <= g_data.max_dist_sq) {
            return true;
        }
    }
    return false; // No witness found for this simplex
}

// Recursively explores all subsets of landmarks to count valid simplices.
void count_simplices_recursive(int* simplex_indices, int current_size, int start_landmark_idx) {
    // A valid simplex must have at least one vertex (dimension 0)
    if (current_size > 0) {
        if (has_witness(simplex_indices, current_size)) {
            g_data.simplex_count++;
        }
    }

    // Stop if the next simplex would exceed the maximum dimension
    if (current_size > g_data.max_simplicial_dimension) {
        return;
    }

    // Explore adding a new landmark to the current simplex
    for (int i = start_landmark_idx; i < g_data.num_landmark_points; ++i) {
        simplex_indices[current_size] = i;
        count_simplices_recursive(simplex_indices, current_size + 1, i + 1);
    }
}


void run_computation() {
    // A k-dimensional simplex has k+1 vertices. Allocate space for the largest possible simplex.
    int* current_simplex_indices = (int*)malloc((g_data.max_simplicial_dimension + 1) * sizeof(int));
    if (!current_simplex_indices) {
        perror("Failed to allocate memory for simplex indices");
        exit(1);
    }

    // Start the recursive process of building and counting simplices.
    // We start with an empty simplex (size 0) from landmark index 0.
    count_simplices_recursive(current_simplex_indices, 0, 0);

    free(current_simplex_indices);
}

void cleanup() {
    // Free landmark points
    for (int i = 0; i < g_data.num_landmark_points; ++i) {
        free(g_data.landmarks[i]);
    }
    free(g_data.landmarks);

    // Free witness points
    for (int i = 0; i < g_data.num_witness_points; ++i) {
        free(g_data.witnesses[i]);
    }
    free(g_data.witnesses);
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
    printf("%llu\n", g_data.simplex_count);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
