#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <string.h>

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

// Benchmark parameters
static int num_points;
static int num_neighbors;
static int target_dimension;
static int ambient_dimension;

// Data and workspace pointers
static double *points;                       // Input point cloud: num_points x ambient_dimension
static double *result_accumulator;           // Stores one value per point to be summed at the end

// Workspace buffers allocated once to be used in run_computation
static double *local_covariance;             // ambient_dimension x ambient_dimension
static double *deflated_matrix;              // ambient_dimension x ambient_dimension
static double *eigenvector;                  // ambient_dimension
static double *temp_vector;                  // ambient_dimension
static double *centroid;                     // ambient_dimension
static int *neighbor_indices;                // num_neighbors
static double *neighbor_dists_sq;            // num_neighbors

// Benchmark result
static double final_result;

// Helper for random double in [0, 1]
double rand_double() {
    return (double)mt_rand() / UINT32_MAX;
}

// Calculates squared Euclidean distance
double distance_sq(const double *p1, const double *p2, int dim) {
    double dist = 0.0;
    for (int i = 0; i < dim; ++i) {
        double diff = p1[i] - p2[i];
        dist += diff * diff;
    }
    return dist;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s num_points num_neighbors target_dim ambient_dim seed\n", argv[0]);
        exit(1);
    }
    num_points = atoi(argv[1]);
    num_neighbors = atoi(argv[2]);
    target_dimension = atoi(argv[3]);
    ambient_dimension = atoi(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    mt_seed(seed);

    // Allocate main data structures
    points = (double *)malloc(num_points * ambient_dimension * sizeof(double));
    result_accumulator = (double *)malloc(num_points * sizeof(double));

    // Allocate workspace buffers
    local_covariance = (double *)malloc(ambient_dimension * ambient_dimension * sizeof(double));
    deflated_matrix = (double *)malloc(ambient_dimension * ambient_dimension * sizeof(double));
    eigenvector = (double *)malloc(ambient_dimension * sizeof(double));
    temp_vector = (double *)malloc(ambient_dimension * sizeof(double));
    centroid = (double *)malloc(ambient_dimension * sizeof(double));
    neighbor_indices = (int *)malloc(num_neighbors * sizeof(int));
    neighbor_dists_sq = (double *)malloc(num_neighbors * sizeof(double));

    if (!points || !result_accumulator || !local_covariance || !deflated_matrix || !eigenvector || !temp_vector || !centroid || !neighbor_indices || !neighbor_dists_sq) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Generate points on a low-dimensional manifold embedded in a high-dimensional space + noise
    for (int i = 0; i < num_points; ++i) {
        double *p = points + i * ambient_dimension;
        // Create base point on a `target_dimension` plane
        for (int j = 0; j < target_dimension; ++j) {
            p[j] = rand_double() * 10.0; // Sample from a hypercube
        }
        for (int j = target_dimension; j < ambient_dimension; ++j) {
            p[j] = 0; // Pad with zeros
        }
        // Add small ambient noise to all dimensions
        for (int j = 0; j < ambient_dimension; ++j) {
            p[j] += (rand_double() - 0.5) * 0.1;
        }
    }
}

void run_computation() {
    const int power_iterations = 10;
    final_result = 0.0;

    for (int i = 0; i < num_points; ++i) {
        // 1. Find k-Nearest Neighbors for point i
        for (int k = 0; k < num_neighbors; ++k) {
            neighbor_dists_sq[k] = HUGE_VAL;
            neighbor_indices[k] = -1;
        }

        for (int j = 0; j < num_points; ++j) {
            if (i == j) continue;

            double d_sq = distance_sq(points + i * ambient_dimension, points + j * ambient_dimension, ambient_dimension);

            int max_dist_k = 0;
            for (int k = 1; k < num_neighbors; ++k) {
                if (neighbor_dists_sq[k] > neighbor_dists_sq[max_dist_k]) {
                    max_dist_k = k;
                }
            }

            if (d_sq < neighbor_dists_sq[max_dist_k]) {
                neighbor_dists_sq[max_dist_k] = d_sq;
                neighbor_indices[max_dist_k] = j;
            }
        }

        // 2. Compute centroid of the neighborhood
        memset(centroid, 0, ambient_dimension * sizeof(double));
        for (int k = 0; k < num_neighbors; ++k) {
            double *neighbor_p = points + neighbor_indices[k] * ambient_dimension;
            for (int d = 0; d < ambient_dimension; ++d) {
                centroid[d] += neighbor_p[d];
            }
        }
        for (int d = 0; d < ambient_dimension; ++d) {
            centroid[d] /= num_neighbors;
        }

        // 3. Compute local covariance matrix
        memset(local_covariance, 0, ambient_dimension * ambient_dimension * sizeof(double));
        for (int k = 0; k < num_neighbors; ++k) {
            double *neighbor_p = points + neighbor_indices[k] * ambient_dimension;
            for (int r = 0; r < ambient_dimension; ++r) {
                double centered_r = neighbor_p[r] - centroid[r];
                for (int c = 0; c < ambient_dimension; ++c) {
                    double centered_c = neighbor_p[c] - centroid[c];
                    local_covariance[r * ambient_dimension + c] += centered_r * centered_c;
                }
            }
        }
        for (int r = 0; r < ambient_dimension * ambient_dimension; ++r) {
            local_covariance[r] /= (num_neighbors > 1 ? num_neighbors - 1 : 1);
        }

        // 4. Find principal components via Power Iteration with Deflation
        memcpy(deflated_matrix, local_covariance, ambient_dimension * ambient_dimension * sizeof(double));
        
        for (int t = 0; t < target_dimension; ++t) {
            // Initialize random vector
            double norm_b = 0.0;
            for (int d = 0; d < ambient_dimension; ++d) {
                eigenvector[d] = rand_double();
                norm_b += eigenvector[d] * eigenvector[d];
            }
            norm_b = sqrt(norm_b);
            for (int d = 0; d < ambient_dimension; ++d) eigenvector[d] /= norm_b;

            // Power Iteration
            for (int iter = 0; iter < power_iterations; ++iter) {
                memset(temp_vector, 0, ambient_dimension * sizeof(double));
                for (int r = 0; r < ambient_dimension; ++r) {
                    for (int c = 0; c < ambient_dimension; ++c) {
                        temp_vector[r] += deflated_matrix[r * ambient_dimension + c] * eigenvector[c];
                    }
                }
                double norm_v = 0.0;
                for (int d = 0; d < ambient_dimension; ++d) norm_v += temp_vector[d] * temp_vector[d];
                norm_v = sqrt(norm_v);
                for (int d = 0; d < ambient_dimension; ++d) eigenvector[d] = temp_vector[d] / norm_v;
            }

            if (t == 0) {
                result_accumulator[i] = eigenvector[0];
            }

            // Deflate matrix: M' = M - lambda * v * v^T
            double lambda = 0.0;
            for(int d=0; d<ambient_dimension; ++d) {
                double mv = 0.0;
                for(int c=0; c<ambient_dimension; ++c) {
                    mv += deflated_matrix[d * ambient_dimension + c] * eigenvector[c];
                }
                lambda += eigenvector[d] * mv;
            }
            
            for (int r = 0; r < ambient_dimension; ++r) {
                for (int c = 0; c < ambient_dimension; ++c) {
                    deflated_matrix[r * ambient_dimension + c] -= lambda * eigenvector[r] * eigenvector[c];
                }
            }
        }
    }

    // 5. Accumulate final result to prevent dead code elimination
    for (int i = 0; i < num_points; ++i) {
        final_result += result_accumulator[i];
    }
}

void cleanup() {
    free(points);
    free(result_accumulator);
    free(local_covariance);
    free(deflated_matrix);
    free(eigenvector);
    free(temp_vector);
    free(centroid);
    free(neighbor_indices);
    free(neighbor_dists_sq);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%f\n", final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}